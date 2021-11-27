#ifndef BASELINE_H
#define BASELINE_H

#include "Abstract.h"
#include "StreamSummary.h"

#define MAXINTER 10

class Baseline : public Abstract{
public:

    Baseline(uint32_t _MEMORY, double ratio = 0.9, std::string _name = "Baseline"){
        NAME = _name;
        MEMORY = _MEMORY;

        LIGHT_LENGTH = MEMORY * ratio * 8;
        bitmap = new BitMap(LIGHT_LENGTH);

        summary = new StreamSummary<ItemPair, COUNT_TYPE>
                (summary->Memory2Size(MEMORY * (1 - ratio)));
    }

    ~Baseline(){
        delete summary;
        delete bitmap;
    }

    void Insert(const ItemPair& item){
        uint32_t bloomPos = (hash(item.item) + item.time) % LIGHT_LENGTH;

        uint32_t interval, tempPos = bloomPos;
        for(interval = 0;interval < MAXINTER;++interval){
            if(bitmap->Get(tempPos)){
                break;
            }
            tempPos = (tempPos + LIGHT_LENGTH - 1) % LIGHT_LENGTH;
        }

        bitmap->Set(bloomPos);

        if(interval >= MAXINTER)
            return;
        ItemPair temp(interval, item.item);

        if(summary->mp->Lookup(temp))
            summary->Add_Data(temp);
        else{
            if(summary->isFull())
                summary->SS_Replace(temp);
            else
                summary->New_Data(temp);
        }

    }

    HashMap Report(COUNT_TYPE HIT){
        return summary->Report(HIT);
    }

private:
    uint32_t LIGHT_LENGTH;

    StreamSummary<ItemPair, COUNT_TYPE>* summary;
    BitMap* bitmap;
};


#endif
