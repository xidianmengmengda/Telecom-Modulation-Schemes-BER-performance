# Telecom: Various Modulation Schemes and their impacts on BER performance
Effects of various modulation schemes on BER performance and data rates: MATLAB simulation

Author: Justin Choi, Nick Jensen, Reif Virtue

The project consists of three parts:
part 2: Testing BPSK in a wireless channel 
part 3: Testing QPSK in a wireless channel and compare it to BPSK
Part 4: Designing an appropriate M-QAM system given a data rate/Bandwidth requirement. Then performing simulations to observe its BER. 



# Overview
This project investigated the suitability of a various modulation schemes including BPSK, QPSK and 16-QAM. 
The simulation found that to meet a data rate of 3Mbps at a 3.125 MHz, a QPSK modulation scheme is sufficient and will meet the requirement. However, if the higher speed demand from the video transmission arises, it is recommended that the island adapts a 16-QAM or higher order modulation scheme to meet the speed requirement. Nevertheless, it is warned that with the use of higher modulation scheme, the modulating signals will need to be used with higher power to ensure enhanced data quality at the receiving end. 

In a real world digital transmission application, different modulations schemes are used for different applications depending on their suitability and requirements.

It becomes clear that while PSK systems are used in low speed communications such as phone lines and fax modems, the QAM and PAM systems are used in much more faster communications such as digital TV and internet modems, as well as mobile communications and w-fi standards. The main reason why QAM is used more widely in high speed transmission than PSK is mainly because since QAM modulates the amplitude as well as the phase, it is more spectrally efficient (Radio Electronics, n.d.). 

In general, higher order modulations such as 64-QAM or 256-QAM has a strong advantage on the data rate speed. They can transmit messages at much higher data rate. This is an intuitive outcome since more symbols mean that there are larger range of bit combinations the symbol schemes can carry.  Hence when the bandwidth is limited, the use of higher order modulation scheme will enhance the speed of the data transmission. This advantage however comes at a price – when the number of symbols increase, the constellation diagram becomes more packed and eventually this leads to the reduction of minimum distance between two adjacent symbols. 

This translates to the fact that increasing the minimum distance will decrease the BER. Therefore, for a lower order modulation scheme such as BPSK, given that the average energy is remains the same as other systems, the minimum distance between the symbols are much larger and hence the BER will be much lower. This theory is confirmed when observing the BER plots of BPSK (part 2) and 16-QAM (Part 4).
