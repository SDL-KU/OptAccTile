# OptAccTile

The optimization process for hardware accelerators with multiple accelerator tiles (OptAccTile):\
\
Optimization of multi-core accelerator performance based on accurate performance estimation, Sunwoo Kim, Youngho Seo, Sungkyung Park, and Chester Sungchung Park, IEEE Access, Feb., 2022.\
\
This algorithm is developed for the performance estimation of multicore accelerators with direct memory access controllers (DMACs). The algorithm predicts a dynamic communication bandwidth of each DMAC based on the runtime state of DMACs, making it possible to estimate the communication amounts handled by DMACs accurately by taking into account the temporal intervals. This algorithm can estimate the performance of a multicore accelerator accurately, regardless of the system communication bandwidth. In addition, the algorithm is used to explore a design space of accelerator core dimensions, and the resulting optimal core dimension provides performance gains compared to the conventional multicore accelerator and single-core accelerator, respectively. This result was also verified by the hardware implementations on Xilinx ZYNQ. 

Contributors: Sunwoo Kim, Youngho Seo and Prof. Chester Sungchung Park\
If you have any questions, please send an email to sunwkim@konkuk.ac.kr\
Original release: Konkuk University, 2022


# How to Run

Multi-core accelerator optimization for convolutional neural networks (CNNs): \
1. Open "opt_tile.m" file in Matlab\

2. Edit layer parameters.\
 %% parameters\
 mac_units_resource: Number of available MAC (multiplication and accumulation) units\
 bram_resource: Maximum available size of on-chip buffer in accelerator (in terms of the numbers of pixels/weights)\
 
 target: Target bandwidth\
 step: Iteration step size for bandwidths\
 
 C_C: Input channels\
 M_C: Output channels\
 E_C: Height of output feature maps\
 F_C: Width of output feature maps\
 R_C: Height of filter weights\
 S_C: Width of filter weights\
 U_C: Stride size\
 
 %% assign layers\
 CNN_Layer: Layer numbers for layer parameters in '%% parameters' section\
 num_clp: Number of cores in the multi-core accelerator\
 
3. Push "F5" to run the Matlab script\

4. Check the optimization results\
 Performance: variable 'exeCycles' is the execution times for accelerator cores in terms of accelerator clock cycles\
 Core dimensions: variable 'A(1)' include the optimiztion results\
  A(1).assign_layers: Layer numbers assigned to accelerator cores\
  A(1).tm: Tile sizes of output channels for accelerstor cores\
  A(1).tc: Tile sizes of input channels for accelerstor cores\
  A(1).te: Tile sizes of OFM height for accelerstor cores\
  A(1).tf: Tile sizes of OFM width for accelerstor cores\
