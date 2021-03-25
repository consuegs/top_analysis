#ifndef DATATRIGGER2016_HPP__
#define DATATRIGGER2016_HPP__

#include "tree/TreeParticles.hpp"
#include "physics.hpp"

namespace dataTrigger2016
{
   bool DataTriggerSelection2016(std::vector<bool> const &diElectronTriggers, std::vector<bool> const &diMuonTriggers, std::vector<bool> const &electronMuonTriggers,
                        std::vector<bool> const &channel, std::vector<bool> const &PD={}, bool const &is2016H=false);
} // namespace selection

#endif /* DATATRIGGER2016_HPP__ */
