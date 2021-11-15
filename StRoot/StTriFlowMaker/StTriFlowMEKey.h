#ifndef StTriFlowMEKey_h
#define StTriFlowMEKey_h
#include <iostream>
#include <fstream>
#include "StRoot/StPicoEvent/StPicoPhysicalHelix.h"
#include "Rtypes.h"
#include <map>
#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"

using namespace std;

class MEKey
{
    public:
        Int_t Centrality;
        Int_t vertex;
        Int_t EventPlane;
        Int_t Event_depth;
        Int_t charge;

        MEKey() {}
        MEKey(Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) : Centrality(i), vertex(j), EventPlane(k), Event_depth(l), charge(m) {}

        bool operator<(const MEKey &r) const
        { // "r == right"
            if(Centrality == r.Centrality)
            {
                if(vertex == r.vertex)
                {
                    if(EventPlane == r.EventPlane)
                    {
                        if(Event_depth == r.Event_depth)
                        {
                            return charge < r.charge;
                        }
                        else
                        {
                            return Event_depth < r.Event_depth;
                        }
                    }
                    else
                    {
                        return EventPlane < r.EventPlane;
                    }
                }
                else
                {
                    return vertex < r.vertex;
                }
            }
            else
            {
                return Centrality < r.Centrality;
            }
        }

        void print()
        {
            cout << "Centrality = " << MEKey::Centrality << endl;
            cout << "vertex = " << MEKey::vertex << endl;
            cout << "EventPlane = " << MEKey::EventPlane << endl;
            cout << "Event_depth = " << MEKey::Event_depth << endl;
            cout << "charge = " << MEKey::charge << endl;
        }

        ClassDef(MEKey, 1)
};

typedef std::map<MEKey, std::vector<Int_t> > vectorIntMap;
typedef std::map<MEKey, std::vector<Float_t> > vectorFloatMap;
typedef std::map<MEKey, std::vector<StPicoPhysicalHelix> > vectorHelixMap;
typedef std::map<MEKey, std::vector<TLorentzVector> > vectorLorentzMap;
typedef std::map<MEKey, std::vector<TVector3> > vectorStThreeFMap;

#endif
