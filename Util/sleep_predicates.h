#pragma once

// Simple Logical ExprEssion Parser
// Author: Andrey Alekseenko <al42and@gmail.com>
// Written as part of MDis project on March 2012

#include <Util/mystl.h>
#include <set>

namespace sleeparser {

template <typename D>
class Expression {
public:
    // Return current expression value for given data
    virtual bool evaluate(const D &data) = 0;
    // Return current expression string representation
    virtual std::string name() const = 0;
    virtual int nargs() const = 0;
    virtual bool setarg(int narg, Expression<D> *val) = 0;
    virtual int priority() const = 0;
};

template <typename D>
class Predicate : public Expression<D> {
public:
    virtual bool add_param(const std::string &s) = 0;
    virtual bool end_params() = 0;
    virtual int nargs() const { return 0; }
    virtual bool setarg(int narg, Expression<D> *val) { return false; }
    virtual int priority() const { return 999; }
};

template <typename D>
class ConstPredicate : public Predicate<D> {
     virtual bool add_param(const std::string &s) { return false; }
     virtual bool end_params() { return true; }
};

template <typename D>
class Operator : public Expression<D> { // Just for distinction
};

// Special predicate for data ranges
// Must have val(const D &data) specified
template <typename D, typename T>
class ListPredicate : public Predicate<D> {
public:
    virtual bool add_param(const std::string &s) { _range.insert(str2any<T>(s)); return true; }
    virtual bool end_params() { return true; }
    virtual bool evaluate(const D &data) { return _range.find(val(data)) != _range.end(); }
    virtual std::string name() const { 
        std::string ret(this->list_name());
        ret += "[ ";
        typename std::set<T>::const_iterator itt = _range.begin();
        for (itt = _range.begin(); itt != _range.end(); ++itt)
            ret += any2str<T>(*itt) + " ";
        return ret + "]";
    };
    virtual std::string list_name() const = 0;
protected:
    virtual T val(const D &data) = 0;
    std::set<T> _range;

};

// Specialize for reading ranges of ints
template<typename D>
class RangePredicate : public ListPredicate<D,int> {
public:
    RangePredicate() : to(false) {};
    virtual bool end_params() { return !to; }
    virtual bool add_param(const std::string &s) {
        if (lower(s) == "to") {
            to = true;
        } else {
            int t = str2any<int>(s);
            if (to)
                for (int i = last + 1; i <= t; ++i)
                    this->_range.insert(i);
            else
                this->_range.insert(t);
            last = t;
            to = false;
        }
        return true;
        // TODO: Also handle <, > etc
    };
private:
    int last;
    bool to;
};

// Special cases of predicates
template <typename D>
class TRUE: public ConstPredicate<D> {
public:
    virtual std::string name() const { return "TRUE";}
    virtual bool evaluate(const D &data) { return true; }
};

template <typename D>
class FALSE: public ConstPredicate<D> {
public:
    virtual std::string name() const { return "FALSE";}
    virtual bool evaluate(const D &data) { return false; }
};

// Specify Operators!

// BiOps
template <typename D>
class BinOp: public Operator<D> {
public:
    virtual int nargs() const { return 2; }
    virtual std::string name() const { 
        return std::string("( ") + this->arg[0]->name() + " " + this->self_name() + " " + this->arg[1]->name() + " )"; 
    }
    virtual bool evaluate(const D &data) { 
        return this->op(arg[0]->evaluate(data), arg[1]->evaluate(data)); 
    }
    virtual bool setarg(int narg, Expression<D> *val) { 
        if (narg > 1) return false;
        arg[narg] = val;
        return true;
    }
    virtual ~BinOp() { delete arg[0]; delete arg[1]; }
protected:
    Expression<D> *arg[2];
    virtual std::string self_name() const = 0;
    virtual bool op(bool a, bool b) = 0;
};

template <typename D>
class UnOp: public Operator<D> {
public:
    virtual int nargs() const { return 1; }
    virtual std::string name() const {return this->self_name() + "( " + this->arg->name() + " )"; }
    virtual bool evaluate(const D &data) { return this->op(arg->evaluate(data)); }
    virtual bool setarg(int narg, Expression<D> *val) { 
        if (narg > 0) return false;
        arg = val;
        return true;
    }
    virtual ~UnOp() { delete arg; }
protected:
    Expression<D> *arg;
    virtual std::string self_name() const = 0;
    virtual bool op(bool a) = 0;
};

template <typename D>
class AND : public BinOp<D> {
public:
    virtual int priority() const { return 30; }
protected:
    virtual std::string self_name() const { return "AND"; };
    virtual bool op(bool a, bool b) {return a && b;}
};

template <typename D>
class OR : public BinOp<D> {
public:
    virtual int priority() const { return 20; }
protected:
    virtual std::string self_name() const { return "OR"; };
    virtual bool op(bool a, bool b) {return a || b;}
};

template <typename D>
class XOR : public BinOp<D> {
public:
    virtual int priority() const { return 40; }
protected:
    virtual std::string self_name() const { return "XOR"; }
    virtual bool op(bool a, bool b) {return (!!a) ^ (!!b);}
};

template <typename D>
class NOT : public UnOp<D> {
public:
    virtual int priority() const { return 80; }
protected:
    virtual std::string self_name() const { return "NOT"; };
    virtual bool op(bool a) {return !a;}
};


template <typename D>
class PARANTHESIS : public Operator<D> { // For (-type paranthesis
public:
    virtual std::string name() const { return std::string("<(>");}
    virtual bool evaluate(const D &data) { return false; }
    virtual int priority() const { return -999; }
    virtual int nargs() const { return 0; }
    virtual bool setarg(int narg, Expression<D> *val) { return false; }
};

///////////////////////////////////// FACTORY
template <typename D>
class ExpressionStringFactory{
public:
    virtual Expression<D>* create(const std::string &name) const {
        if      (name == "and")
            return new AND<D>;
        else if (name == "or" )
            return new OR<D>;
        else if (name == "xor")
            return new XOR<D>;
        else if (name == "not")
            return new NOT<D>;
        else if (name == "true")
            return new TRUE<D>;
        else if (name == "false")
            return new FALSE<D>;
        else return NULL;
    }
    virtual bool is_operator(const std::string &name) const {
        return (name == "and" || name == "or" || name == "xor" || name == "not" || name == "true" || name == "false");
    }
    virtual bool is_keyword(const std::string &name) const { return (this->is_operator(name) || name == "(" || name == ")"); }
};




}; // namespace sleep
