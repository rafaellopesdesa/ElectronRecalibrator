#ifndef EP_COMBINATION_TOOL_HPP
#define EP_COMBINATION_TOOL_HPP

#include <string>
#include "GBRForest.h"
#include "../../CMS3.h"

class EpCombinationTool

{
    public:
        EpCombinationTool();
        ~EpCombinationTool();
        // forbid copy and assignment, since we have a custom deleter
        EpCombinationTool(const EpCombinationTool &other) = delete;
        EpCombinationTool & operator=(const EpCombinationTool &other) = delete;

        bool init(const CORE_GBR::GBRForest *forest) ;
        bool init(const std::string& regressionFile, const std::string& bdtName);
	void combine(int elsIdx, float old_energy, float old_energy_error, float& recalib_energy) const;


    private:
        const CORE_GBR::GBRForest* m_forest;
        bool  m_ownForest;

};


#endif
