#pragma once

#include <Core/topology.h>
#include <string>

class AtomFilter {
    public:
        AtomFilter(const std::string &expr);
        ~AtomFilter();
        bool operator()(const Atom &a) { return this->process(a); }
        bool process(const Atom &a);
        void masquerade(const Atom *atoms, size_t size, int *out, int mask);
        std::string str() const;
    private:
        // I don't want to include any SLEeP-related headers here, since they are heavily-templated ones, so it would make compilation even slower
        struct _parser_wrapper;
        _parser_wrapper *pw;
};


std::ostream& operator<<(std::ostream& str, AtomFilter const& data);
