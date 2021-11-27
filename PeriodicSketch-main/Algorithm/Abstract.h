#ifndef ABSTRACT_H
#define ABSTRACT_H

#include "hash.h"
#include "Util.h"

class Abstract{
public:
    std::string NAME;
    uint64_t MEMORY;

    Abstract(){}
    virtual ~Abstract(){};

    virtual void Insert(const ItemPair& item) = 0;
    virtual HashMap Report(COUNT_TYPE HIT) = 0;
};

#endif
