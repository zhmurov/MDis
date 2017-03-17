#include "atomfilter.h"
#include <Util/sleep.h>
#include <Util/wrapper.h>
#include <exception>
#include <Core/topology.h>
#include <string>


// Define all predicates
class P_resid: public sleeparser::RangePredicate<Atom> {
    virtual int val(const Atom &data) { return data.resid; }
    virtual std::string list_name() const { return "RESID";}
};
class P_index: public sleeparser::RangePredicate<Atom> {
    virtual int val(const Atom &data) { return data.id; }
    virtual std::string list_name() const { return "INDEX";}
};
class P_name: public sleeparser::ListPredicate<Atom, std::string> {
    virtual std::string val(const Atom &data) { return std::string(data.name); }
    virtual std::string list_name() const { return "NAME";}
};
class P_chain: public sleeparser::ListPredicate<Atom, char> {
    virtual char val(const Atom &data) { return data.chain; }
    virtual std::string list_name() const { return "CHAIN";}
};
class P_segid: public sleeparser::ListPredicate<Atom, std::string> {
    virtual std::string val(const Atom &data) { return std::string(data.segment); }
    virtual std::string list_name() const { return "SEGID";}
};


class AtomPredicateFactory : public sleeparser::ExpressionStringFactory<Atom> {
    virtual sleeparser::Expression<Atom>* create(const std::string &name) const {
        if      (name == "resid") return new P_resid;
        else if (name == "index") return new P_index;
        else if (name == "name" ) return new P_name ;
        //else if (name == "chain") return new P_chain; // Too errorprone!
        else if (name == "segid") return new P_segid;
        else return sleeparser::ExpressionStringFactory<Atom>().create(name);
    }
};

struct AtomFilter::_parser_wrapper {
    sleeparser::Parser<Atom> *parser;
};

AtomFilter::AtomFilter(const std::string &expr) : pw(new _parser_wrapper()) {
    try {
        this->pw->parser = new sleeparser::Parser<Atom>(expr, AtomPredicateFactory());
    } catch (const char *e) {
        DIE("Error parsing '%s': %s", expr.c_str(), e);
    }
}

AtomFilter::~AtomFilter() { 
    delete this->pw->parser;
    delete this->pw;
}

bool AtomFilter::process(const Atom &a) { 
    try {
        return (*pw->parser)(a); 
    } catch (const char *e) {
        DIE("Error processing atom: %s", e);
    }
}

void AtomFilter::masquerade(const Atom *atoms, size_t size, int *out, int mask) {
    for (size_t i = 0; i < size; ++i) {
        if (this->process(atoms[i]))
            out[i] = mask;
    }
}

std::string AtomFilter::str() const { return pw->parser->str(); }

 
std::ostream& operator<<(std::ostream& str, AtomFilter const& data) {
    str << data.str();
    return str;
}

