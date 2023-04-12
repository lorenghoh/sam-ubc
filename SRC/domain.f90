! Set the domain dimensionality, size and number of subdomains.

module domain

       integer, parameter :: YES3D = 1  ! Domain dimensionality: 1 - 3D, 0 - 2D
       integer, parameter :: nx_gl = 64 ! Number of grid points in X
       integer, parameter :: ny_gl = 64 ! Number of grid points in Y
       integer, parameter :: nz_gl = 64 ! Number of pressure (scalar) levels
       integer, parameter :: nsubdomains_x  = 4 ! No of subdomains in x
       integer, parameter :: nsubdomains_y  = 8 ! No of subdomains in y


       ! define # of points in x and y direction to average for 
       !   output relating to statistical moments.
       ! For example, navgmom_x = 8 means the output will be   
       !  8 times coarser grid than the original.
       ! If don't wanna such output, just set them to -1 in both directions. 
       ! See Changes_log/README.UUmods for more details.
       integer, parameter :: navgmom_x = -1 
       integer, parameter :: navgmom_y = -1 

       integer, parameter :: ntracers = 0 ! number of transported tracers (dotracers=.true.)
       
! Note:
!  * nx_gl and ny_gl should be a factor of 2,3, or 5 (see User's Guide)
!  * if 2D case, ny_gl = nsubdomains_y = 1 ;
!  * nsubdomains_x*nsubdomains_y = total number of processors
!  * if one processor is used, than  nsubdomains_x = nsubdomains_y = 1;
!  * if ntracers is > 0, don't forget to set dotracers to .true. in namelist 

end module domain


! Acceptable nx_gl or ny_gl for the FFT to work
!
!           4           5           6           8           9          10          12          15
!          15          16          18          20          24          25          27          30
!          30          32          36          40          45          48          50          54
!          54          60          64          72          75          80          81          90
!          90          96         100         108         120         125         128         135
!         135         144         150         160         162         180         192         200
!         200         216         225         240         243         250         256         270
!         270         288         300         320         324         360         375         384
!         384         400         405         432         450         480         486         500
!         500         512         540         576         600         625         640         648
!         648         675         720         729         750         768         800         810
!         810         864         900         960         972        1000        1024        1080
!        1080        1125        1152        1200        1215        1250        1280        1296
!        1296        1350        1440        1458        1500        1536        1600        1620
!        1620        1728        1800        1875        1920        1944        2000        2025
!        2025        2048        2160        2187        2250        2304        2400        2430
!        2430        2500        2560        2592        2700        2880        2916        3000
!        3000        3072        3125        3200        3240        3375        3456        3600
!        3600        3645        3750        3840        3888        4000        4050        4096
!        4096        4320        4374        4500        4608        4800        4860        5000
!        5000        5120        5184        5400        5625        5760        5832        6000
!        6000        6075        6144        6250        6400        6480        6561        6750
!        6750        6912        7200        7290        7500        7680        7776        8000
!        8000        8100        8192        8640        8748        9000        9216        9375
!        9375        9600        9720       10000

