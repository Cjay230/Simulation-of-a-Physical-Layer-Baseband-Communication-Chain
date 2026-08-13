This project implements and evaluates, in MATLAB, the complete transmitter–receiver chain of a physical-layer 
baseband digital communication system. An analog message m(t) was converted into a digital bit stream through 
sampling, quantization, and source coding (Part I), transmitted over a noisy channel using digital modulation, and 
recovered through demodulation, detection, and the corresponding decoding stages (Part II). Each block was 
implemented as a MATLAB local function with the required signature, verified in isolation. 
Four modules were carried out in order: (1) signal representation through sampling and sinc reconstruction, including 
a demonstration of aliasing below the Nyquist rate; (2) two-level and uniform multi-level quantization with mean
square-error (MSE) analysis for M = 4, 16, and 32; (3) source coding of the quantized symbol stream, comparing a 
fixed-length baseline against an optimal Huffman code built from the empirical symbol distribution and benchmarked 
against the source entropy; and (4) BPSK and QPSK modulation over an additive white Gaussian noise (AWGN) 
channel, with bit-error-rate (BER) estimation over 10⁶ bits and comparison to the theoretical Q-function curve.
