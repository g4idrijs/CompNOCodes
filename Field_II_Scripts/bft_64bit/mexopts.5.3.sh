#
# mexopts.sh   Shell script for configuring MEX-file creation script,
#               mex.
#
# usage:        Do not call this file directly; it is sourced by the
#               mex shell script.  Modify only if you don't like the
#               defaults after running mex.  No spaces are allowed
#               around the '=' in the variable assignment.
#
# Copyright (c) 1992-95 by The MathWorks, Inc.
# $Revision: 1.23 $  $Date: 1996/10/24 18:08:03 $
#----------------------------------------------------------------------------
#
    case "$Arch" in
        Undetermined)
#----------------------------------------------------------------------------
# Change this line if you need to specify the location of the MATLAB
# root directory.  The cmex script needs to know where to find utility
# routines so that it can determine the architecture; therefore, this
# assignment needs to be done while the architecture is still
# undetermined.
#----------------------------------------------------------------------------
            MATLAB="$MATLAB"
            ;;
        alpha)
#----------------------------------------------------------------------------
            CC='cc'
            CFLAGS='-ieee_with_no_inexact -std1 -assume noaccuracy_sensitive'
            CLIBS=''
            COPTIMFLAGS='-O'
            CDEBUGFLAGS='-g'
#
            FC='f77'
            FFLAGS='-shared'
            FLIBS='-lUfor -lfor -lFutil'
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD='ld'
            LDFLAGS="-expect_unresolved '*' -shared -hidden -exported_symbol $ENTRYPOINT -exported_symbol mexVersion"
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
        hp700)
#----------------------------------------------------------------------------
            CC='CC'
            CFLAGS='+z -D_HPUX_SOURCE +a1'
            CLIBS=''
            COPTIMFLAGS='-O'
            CDEBUGFLAGS='-g'
#
            FC='f77'
            FFLAGS='+z'
            FLIBS=''
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD='ld'
            LDFLAGS="-b +e $ENTRYPOINT +e mexVersion"
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
        ibm_rs)
#----------------------------------------------------------------------------
            CC='cc'
            CFLAGS='-qlanglvl=ansi'
            CLIBS='-lm'
            COPTIMFLAGS='-O'
            CDEBUGFLAGS='-g'
#
            FC='f77'
            FFLAGS=''
            FLIBS="$MATLAB/extern/lib/ibm_rs/fmex1.o -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD='cc'
            LDFLAGS="-bI:$MATLAB/extern/lib/ibm_rs/exp.ibm_rs -bM:SRE -e $ENTRYPOINT"
            LDOPTIMFLAGS='-s'
            LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
        lnx86)
#----------------------------------------------------------------------------
            CC='gcc'
            CFLAGS='-DCPU=686 -D__AUTO_DISABLE -Wall'
            CLIBS='-L/usr/i486-linux-libc5/lib'
            COPTIMFLAGS='-O3'
            CDEBUGFLAGS='-g'
#
# Use these flags for using f2c and gcc for Fortan MEX-Files
#
            FC='f2c'
            FOPTIMFLAGS=''
            FFLAGS=''
            FDEBUGFLAGS='-g'
            FLIBS=''
#
# Use these flags for using the Absoft F77 Fortran Compiler
#
        #   FC='f77'
        #   FOPTIMFLAGS=''
        #   FFLAGS='-f -N9 -N70'
        #   FDEBUGFLAGS='-gg'
        #   FLIBS='-lf77'
#
            LD='gcc'
            LDFLAGS='-shared -rdynamic'
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
        sgi)
#----------------------------------------------------------------------------
            CC='cc'
            CFLAGS='-ansi'
            CLIBS=''
            COPTIMFLAGS='-O'
            CDEBUGFLAGS='-g'
#
            FC='f77'
            FFLAGS=''
            FLIBS=''
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD='ld'
            LDFLAGS="-shared -U -Bsymbolic -exported_symbol $ENTRYPOINT -exported_symbol mexVersion"
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''
            ;;
#----------------------------------------------------------------------------
        sgi64)
# R8000 only: The default action of mex is to generate full MIPS IV
#             (R8000) instruction set.
#----------------------------------------------------------------------------
            CC='cc'
            CFLAGS='-mips4'
            CLIBS=''
            COPTIMFLAGS='-O'
            CDEBUGFLAGS='-g'
#
            FC='f77'
            FFLAGS=''
            FLIBS=''
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD='ld'
            LDFLAGS="-mips4 -shared -U -Bsymbolic -exported_symbol $ENTRYPOINT -exported_symbol mexVersion"
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''
            ;;
#----------------------------------------------------------------------------
        sol2)
#----------------------------------------------------------------------------
            CC='cc'
            CFLAGS=''
            CLIBS=''
            COPTIMFLAGS='-O'
            CDEBUGFLAGS='-g'
#
            FC='f77'
            FFLAGS='-G'
            FLIBS=''
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD='ld'
            LDFLAGS="-G -M $MATLAB/extern/lib/sol2/$MAPFILE"
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
        sun4)
#----------------------------------------------------------------------------
# A dry run of the appropriate compiler is done in the mex script to
# generate the correct library list. Use -v option to see what
# libraries are actually being linked in.
#----------------------------------------------------------------------------
            CC='acc'
            CFLAGS='-DMEXSUN4'
            CLIBS="$MATLAB/extern/lib/sun4/libmex.a -lm"
            COPTIMFLAGS='-O'
            CDEBUGFLAGS='-g'
#
            FC='f77'
            FFLAGS=''
            FLIBS='$MATLAB/extern/lib/sun4/libmex.a -lm'
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD='ld'
            LDFLAGS='-d -r -u _mex_entry_pt -u _mexFunction'
            LDOPTIMFLAGS='-x'
            LDDEBUGFLAGS=''
#----------------------------------------------------------------------------
            ;;
    esac
#############################################################################
#
# Architecture independent lines:
#
#     Set and uncomment any lines which will apply to all architectures.
#
#----------------------------------------------------------------------------
#           CC="$CC"
#           CFLAGS="$CFLAGS"
#           COPTIMFLAGS="$COPTIMFLAGS"
#           CDEBUGFLAGS="$CDEBUGFLAGS"
#           CLIBS="$CLIBS"
#
#           FC="$FC"
#           FFLAGS="$FFLAGS"
#           FOPTIMFLAGS="$FOPTIMFLAGS"
#           FDEBUGFLAGS="$FDEBUGFLAGS"
#           FLIBS="$FLIBS"
#
#           LD="$LD"
#           LDFLAGS="$LDFLAGS"
#           LDOPTIMFLAGS="$LDOPTIMFLAGS"
#           LDDEBUGFLAGS="$LDDEBUGFLAGS"
#----------------------------------------------------------------------------
#############################################################################
