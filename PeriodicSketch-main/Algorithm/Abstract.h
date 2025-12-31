#ifndef ABSTRACT_H
#define ABSTRACT_H

#include "../common/hash.h"
#include "../common/Util.h"

class Abstract{
public:
    std::string NAME;
    uint64_t MEMORY;
    double   RATIO;

    Abstract(){}
    virtual ~Abstract(){};

    virtual void Insert(const ItemPair& item) = 0;
    //virtual HashMap Report(COUNT_TYPE HIT) = 0;
    virtual HashMap Report1(COUNT_TYPE HIT) = 0;
};

#endif
