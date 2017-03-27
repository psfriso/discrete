!
! CommLine & I/O
!
   integer, parameter :: NFILES = 12
 
   integer unit_i, unit_o, unit_top, unit_r, unit_ener, unit_traj, unit_pdb, unit_v, unit_l
 
   type(commLineOption) :: files(NFILES) = (/&
   commLineOption("-i",    "param",          "formatted",   "old",     "Settings"),&
   commLineOption("-top",  "topology",       "unformatted", "old",     "Topology"),&
   commLineOption("-r",    "coordinates",    "unformatted", "old",     "Initial coordinates"),&
   commLineOption("-ener", "energy",         "formatted",   "unknown", "Energies"),&
   commLineOption("-traj", "trajectory.pdb", "formatted",   "unknown", "Trajectory (PDB)"),&
   commLineOption("-o",    "log",            "formatted",   "unknown", "Calculation Log"),&
   commLineOption("-x",    "trajectory.crd", "formatted",   "unknown", "Trajectory (CRD)"),&   
   commLineOption("-rst",  "restart",        "unformatted", "unknown", "Restart coordinates"), &
   commLineOption("-rstv", "restartVel",     "unformatted", "unknown", "Restart velocities"), &
   commLineOption("-v",    "velocities",     "unformatted", "old",     "Initial velocities"), &
   commLineOption("-sdf",  "sysdef",         "formatted",   "old",     "System Definition"), &
   commLineOption("-l",    "log",            "formatted",   "unknown",  "Calculation log, machine readable") &
   /)        
!
   type(systemDef) :: sysDef
   type(strucTopology), allocatable, target :: refTop(:), simTop(:)
   type(simulRef), allocatable :: refStruc(:)
   type(simulRT), allocatable :: simStruc(:)
!   
   real calcEkin
   integer ibloc, irepl
   real*8 tacum, tinit, tsetup, tfin
!
   integer i,j,ns