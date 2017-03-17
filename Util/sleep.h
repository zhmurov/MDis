#pragma once

// Simple Logical ExprEssion Parser
// Author: Andrey Alekseenko <al42and@gmail.com>
// Written as part of MDis project on March 2012
//
//
//
// Bicycle bicycle bicycle
// I want to ride my bicycle bicycle bicycle
//
// I want to ride my bicycle
// I want to ride my bike
// I want to ride my bicycle
// I want to ride it where I like
//

#include <map>
#include <stack>
#include <queue>
#include <string>
#include <sstream>
#include <iostream>
#include "sleep_predicates.h"

namespace sleeparser {

template <typename D>
void push_expr(std::stack< Expression<D>* > &stack, Expression<D> *op) {
    const int nargs = op->nargs();
    for (int i = 0; i < nargs; ++i) {
        if (!op->setarg(i, stack.top()))
            throw("Unable to set argument");
        stack.pop();
    }
    stack.push(op);
}

template <typename D>
class Parser {
public:
    Parser(const std::string &expr, const ExpressionStringFactory<D> &esf ) {
        // Adding spaces arount braces to simplify tokenization
        std::string newexpr = string_replace(expr, "(", " ( ");
        newexpr = string_replace(newexpr, ")", " ) ");
        newexpr = std::string("( ") + newexpr + " )";

        // Tokenize!
        std::istringstream iss(newexpr);
        std::string token;
        std::queue<std::string> tokens;
        while (iss >> token)
            tokens.push(token);

        // Do Shunting-yard'ing
        std::stack<Operator<D>*> stack;
        std::stack<Expression<D>*> output;
        while( ! tokens.empty() ) {
            token = tokens.front();
            std::string ltoken = lower(token);
            if (esf.is_operator(ltoken)) {
                Operator<D> *t;
                t = dynamic_cast< Operator<D>* > ( esf.create(ltoken) );
                while ( !stack.empty() ) {
                    if (stack.top()->priority() >= t->priority()) {
                        push_expr(output, stack.top());
                        stack.pop();
                    } else break;
                }
                stack.push(t);
                tokens.pop();
            } else if (ltoken == "(") {
                stack.push( new PARANTHESIS<D> );
                tokens.pop();
            } else if (ltoken == ")") {
                Expression<D> *t;
                bool ok = false;
                while (!stack.empty()) {
                    t = stack.top();
                    if (t->priority() != PARANTHESIS<D>().priority()) { // We consider same paranthesis types have same priority
                        push_expr(output, t);
                    } else {
                        ok = true;
                        break;
                    }
                    stack.pop();
                }
                if (!ok) throw "Unbalanced parantheses!";
                tokens.pop();
            } else { // Custom predicate
                Predicate<D> *t;
                t = dynamic_cast< Predicate<D>* > ( esf.create(ltoken) );
                if (t == NULL) // create custom function
                    throw "Unknown predicate";
                while( true ) {
                    tokens.pop();
                    if (tokens.empty()) break;
                    token = tokens.front();
                    ltoken = lower(token);
                    if (!esf.is_keyword(ltoken)) {
                        if (!t->add_param(token))
                            throw "Can not add parameter to predicate";
                    } else
                        break;
                }
                if (!t->end_params())
                    throw "Can not stop adding parameters to predicate";
                push_expr(output, t);
            }
        }
        if (output.size() != 1)
            throw "Wrong output size";
        // Build parsing tree
        _root = output.top();
    }
    bool operator()(const D &data) {
        return _root->evaluate(data);
    }
    std::string str() const {
        return _root->name();
    }
private:
    Expression<D> *_root;
};

}; // namespace sleep

template <typename D>
std::ostream& operator<<(std::ostream& str, sleeparser::Parser<D> const& data)
{
        str << "{ " << data.str() << " }";
        return str;
}


