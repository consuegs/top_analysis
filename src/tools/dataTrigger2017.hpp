#ifndef DATATRIGGER2017_HPP__
#define DATATRIGGER2017_HPP__

#include "tree/TreeParticles.hpp"
#include "physics.hpp"

namespace dataTrigger2017
{
   bool DataTriggerSelection2017(std::vector<bool> const &diElectronTriggers, std::vector<bool> const &diMuonTriggers, std::vector<bool> const &electronMuonTriggers,
                        std::vector<bool> const &channel, std::vector<bool> const &PD={}, bool const &is2017AB=false);
} // namespace selection

#endif /* DATATRIGGER2017_HPP__ */
