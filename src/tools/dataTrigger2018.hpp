#ifndef DATATRIGGER2018_HPP__
#define DATATRIGGER2018_HPP__

#include "tree/TreeParticles.hpp"
#include "physics.hpp"

namespace dataTrigger2018
{
   bool DataTriggerSelection2018(std::vector<bool> const &diElectronTriggers, std::vector<bool> const &diMuonTriggers, std::vector<bool> const &electronMuonTriggers,
                        std::vector<bool> const &channel, std::vector<bool> const &PD={});
} // namespace selection

#endif /* DATATRIGGER2018_HPP__ */
