# OptAccTile

The optimization process for hardware accelerators with multiple accelerator tiles ("OptAccTile):\
\
Optimization of multi-core accelerator performance based on accurate performance estimation, Sunwoo Kim, Youngho Seo, Sungkyung Park, and Chester Sungchung Park, IEEE Access, Feb., 2022.\
\
This algorithm is developed for the performance estimation of multicore accelerators with direct memory access controllers (DMACs). The algorithm predicts a dynamic communication bandwidth of each DMAC based on the runtime state of DMACs, making it possible to estimate the communication amounts handled by DMACs accurately by taking into account the temporal intervals. This algorithm can estimate the performance of a multicore accelerator accurately, regardless of the system communication bandwidth. In addition, the algorithm is used to explore a design space of accelerator core dimensions, and the resulting optimal core dimension provides performance gains compared to the conventional multicore accelerator and single-core accelerator, respectively. This result was also verified by the hardware implementations on Xilinx ZYNQ. 

Contributors: Sunwoo Kim, Youngho Seo and Prof. Chester Sungchung Park\
If you have any questions, please send an email to sunwkim@konkuk.ac.kr\
Original release: Konkuk University, 2022\
