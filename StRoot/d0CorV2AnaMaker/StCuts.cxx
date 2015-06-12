/* **************************************************
 *  Cuts namespace.
 *
 *  Authors:  Xin Dong        (xdong@lbl.gov)
 *            Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *            Jochen Thaeder  (jmthader@lbl.gov)   
 *
 * **************************************************
 */

#include "StCuts.h"

namespace mycuts
{
   // path to lists of triggers prescales
   // lists are obtained from http://www.star.bnl.gov/protected/common/common2014/trigger2014/plots_au200gev/
   std::string const prescalesFilesDirectoryName = "./run14AuAu200GeVPrescales";
   //event
   UShort_t const triggerWord = 0x1F; //first five bits see http://rnc.lbl.gov/~xdong/SoftHadron/picoDst.html
   float const vz = 6.0;// cm.
   float const vzVpdVz = 3.0; // 3 cm.

   //tracking
   int const nHitsFit = 20;
   bool const requireHFT = true;
   float const minPt = 0.6;

   //pions
   float const nSigmaPion = 3.0;

   //kaons
   float const nSigmaKaon = 2.0;
   float const kTofBetaDiff = 0.03;

   // tree kaonPion pair cuts
   float const cosTheta = 0.995; // minimum
   float const dcaDaughters = 0.0050; // maximum
   float const decayLength = 0.0030; // minimum
   float const minMass = 1.6;
   float const maxMass = 2.1;

   // histograms kaonPion pair cuts
   float const qaNHitsFit = 20;
   float const qaNSigmaKaon = 2.0;
   float const kDca = 0.008; // minimum
   float const pDca = 0.008;
   // Hadron cuts
   float const hadronPtMin= 0.2;
   float const hadronPtMax = 2.0;
   float const corDetaMin = 0.75;
   float const corDetaMax = 10;
}
