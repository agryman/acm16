# Coded by TAW, 2008-2012.
# Version 16 (06 Apr 2012). [Enabling choice of lambda switching,
#                            and various basis choices,
#                            Extra zero AM ops via Racah coeffs.]

# This is an attempt to calculate matrix elements for ACM stuff,
# with the factorisation of states into a direct product of
# beta states (with parameter lambda) labelled by (nu),
# and SO(5) spherical harmonics, labelled by (v,alpha,L,M).
# We ignore the M throughout by dealing with reduced matrix elements.

# This code is being honed alongside an upcoming publication by
# T.A. Welsh and D.J. Rowe.


# It was based on previous publications:

# A pretty full explanation in Chapter 4 of the book [RowanWood]
#   "Fundamentals of Nuclear Models: Foundational Models"
#   by D.J. Rowe and J.L. Wood,
#   World Scientific (Singapore), 2010.

# Also see the paper [RWC2009]
#   "The Bohr model a
#   by D.J. Rowe, T.A. Welsh and M.A. Caprio,
#   Phys. Rev. C79 (2009) 054304. 

# An important precursor to this is the paper [Rowe2004]
#   "A computationally tractable version of the collective model"
#   by D.J. Rowe, Nucl. Phys. A 735 (2004) 372-392.


# At the end, the code below also contains a couple of routines that
# aid the production of the data for the particular Hamiltonians
# considered in [RWC2009].

###########################################################################
###########################################################################

# We make extensive use of the LinearAlgebra library. In particular,
# this provides the diagonalisation procedure that we use.

with(LinearAlgebra):

###########################################################################
####----------------------- Global Constants --------------------------####
###########################################################################

# Here we specify various constants.
# They should not be altered by the user.
# They are mainly used to signify certain operators.
# Hamiltonians and other operators will be expressed in terms of
# these values.

# The following ten operators are the "basic" radial operators.
# The way that they alter lambda is not fixed, but is determined
# automatically. In particular, they will eventually be exchanged
# for operators in which the shift is specific.
# In particular, the final three (which must be numbered consecutively)
# will be replaced by other operators which specify a shift of -1,0 or +1.
# The quadratic operators may (or not) also be exchanged for such.

Radial_Sm:=-1:          #  SU(1,1) operator S-
Radial_S0:=0:           #  SU(1,1) operator S0
Radial_Sp:=1:           #  SU(1,1) operator S+
Radial_b2:=2:           #  beta^2
Radial_bm2:=3:          #  1/beta^2
Radial_D2b:=4:          #  d^2/d(beta)^2
Radial_bDb:=5:          #  beta*d/d(beta)
Radial_b:=6:            # beta
Radial_bm:=7:           # 1/beta
Radial_Db:=8:           # d/d(beta)

# the list of basic radial operators (the order here is important later).

Radial_List:=[Radial_Sm, Radial_S0, Radial_Sp, Radial_b2, Radial_bm2,
              Radial_D2b, Radial_bDb, Radial_b, Radial_bm, Radial_Db]:


# The next set of nine operators are not basic in that they stipulate
# how lambda is to change under their application. They each pertain
# to one of the three odd degree operators beta, 1/beta, d/d(beta).

# the next three, which should be consecutive, shift lambda by -1.

Radial_b_ml:=21:        # beta (with lambda--)
Radial_bm_ml:=22:       # 1/beta (with lambda--)
Radial_Db_ml:=23:       # d/d(beta) (with lambda--)

# the next three, which should be consecutive, shift lambda by 0.
# These will be obtained analytically, involving the taking the
# square root of a degree 2 operator matrix.

Radial_b_sq:=24:        # beta obtained through sqrt(beta^2),
Radial_bm_sq:=25:       # beta^(-1) obtained through sqrt(beta^(-2)),
                        # (both with no shifting of lambda)
Radial_Db_sq:=26:       # d/d(beta) obtained through sqrt(beta^(-2)),

# the next three, which should be consecutive, shift lambda by +1.

Radial_b_pl:=27:        # beta (with lambda++)
Radial_bm_pl:=28:       # 1/beta (with lambda++)
Radial_Db_pl:=29:       # d/d(beta) (with lambda++)

# Differences between above sets

Radial_shift_ml:=Radial_b_ml-Radial_b:
Radial_shift_sq:=Radial_b_sq-Radial_b:
Radial_shift_pl:=Radial_b_pl-Radial_b:

# The following requires an integer parameter.
# This parameter is given as the next item in the list, which is
# itself followed by Radial_flag (useful for reading from right).

Radial_id:=50:      # identity operator for shifting lambda by 2*k

Radial_flag:=-10:   # signals parameter for above operator


# Now for operators on the spherical space: these are all the spherical
# harmonics for v<=6.
# There should be SO5 CG coefficients files present for each of these.

SpHarm_112:=1102;       # denotes Y_{112}
SpHarm_212:=2102;       # denotes Y_{212}
SpHarm_214:=2104;       # denotes Y_{214}
SpHarm_310:=3100;       # denotes Y_{310}
SpHarm_313:=3103;       # denotes Y_{313}
SpHarm_314:=3104;       # denotes Y_{314}
SpHarm_316:=3106;       # denotes Y_{316}
SpHarm_412:=4102; SpHarm_414:=4104; SpHarm_415:=4105;
SpHarm_416:=4106; SpHarm_418:=4108;
SpHarm_512:=5102; SpHarm_514:=5104; SpHarm_515:=5105; SpHarm_516:=5106;
SpHarm_517:=5107; SpHarm_518:=5108; SpHarm_51A:=5110;
SpHarm_610:=6100; SpHarm_613:=6103; SpHarm_614:=6104; SpHarm_616:=6106;
SpHarm_626:=6206; SpHarm_617:=6107; SpHarm_618:=6108; SpHarm_619:=6109;
SpHarm_61A:=6110; SpHarm_61C:=6112; 

# Use the following table to access the corresponding [v,alpha,L].
# (could use iquo(#,1000) etc.)

SpHarm_Table:=table([
  SpHarm_112=[1,1,2],
  SpHarm_212=[2,1,2], SpHarm_214=[2,1,4],
  SpHarm_310=[3,1,0], SpHarm_313=[3,1,3], SpHarm_314=[3,1,4],
  SpHarm_316=[3,1,6],
  SpHarm_412=[4,1,2], SpHarm_414=[4,1,4], SpHarm_415=[4,1,5],
  SpHarm_416=[4,1,6], SpHarm_418=[4,1,8],
  SpHarm_512=[5,1,2], SpHarm_514=[5,1,4], SpHarm_515=[5,1,5],
  SpHarm_516=[5,1,6], SpHarm_517=[5,1,7], SpHarm_518=[5,1,8],
  SpHarm_51A=[5,1,10],
  SpHarm_610=[6,1,0], SpHarm_613=[6,1,3], SpHarm_614=[6,1,4],
  SpHarm_616=[6,1,6], SpHarm_626=[6,2,6], SpHarm_617=[6,1,7],
  SpHarm_618=[6,1,8], SpHarm_619=[6,1,9], SpHarm_61A=[6,1,10],
  SpHarm_61C=[6,1,12] 
]):

# Form a list of the available indices in this table.

SpHarm_Indices:=map(op,[indices(SpHarm_Table)]):


# The following indicates which Spherical Harmonics are currently
# available when constructing representation matrices on the full
# SU(1,1) x SO(5) space using RepXspace_Prod().
# An error results if other SpHarms are used in that procedure.
# [There is no reason why this shouldn't be extended to the set above,
#  apart from perhaps having to put them all in the writeup!]
# However, for representations purely on the SO5 space, this
# is not used. If the matrix elements are not then available, an
# error will result for another reason (possibly file doesn't exist).

SpHarm_List:=[SpHarm_112,SpHarm_212,SpHarm_214,
              SpHarm_310,SpHarm_313,SpHarm_314,SpHarm_316]:

# The following operators intrinsically affect the whole product space

Xspace_Pi:=100:      # For operator  pi;
Xspace_PiPi2:=101:   # For operator  [pi x pi]_{v=2,L=2};
Xspace_PiPi4:=102:   # For operator  [pi x pi]_{v=2,L=4};
Xspace_PiqPi:=120:   # For the "grubby" operator  [q x pi x pi]_{v=3,L=0}.

# The following two only affect the spherical space.

SpDiag_sqLdim:=900:  # Provides a diagonal matrix with entries
                     # (-1)^{L_i}*sqrt(2L_i+1)
SpDiag_sqLdiv:=901:  # Provides a diagonal matrix with entries
                     # (-1)^{L_i}/sqrt(2L_i+1)

# the full list of spherical operators..

Spherical_List:=[op(SpHarm_Indices),SpDiag_sqLdim,SpDiag_sqLdiv]:

# Use the following to determine when operators act solely on
# the radial space (<=Radial_Max), solely on the spherical space
# (>=Spherical_Min), or act intrinsically on the cross-product space.
# Operators in these ranges need not actually exist.

Radial_Max:=99:
Spherical_Min:=900:



# The following provide conversion factors from the spherical harmonics
# to more physically relevant operators
# (Note that often (e.g. by RepSO5_Y), the operator will be represented
# with the 4*Pi already incorporated - and the FourPi should be cancelled).
# Note that evalf will need to be used somewhere further down the line.

FourPi:=4*Pi;
Convert_red:=1/FourPi;         # converts ME_SO5red to <v3|||v2|||v1>

Convert_112:=FourPi/sqrt(15);       # multiplies Y112 to get Q
Convert_212:=-FourPi*sqrt(2/105);   # multiplies Y212 to get [QxQ]_(L=4)
Convert_310:=FourPi/3;              # multiplies Y310 to get cos(3*gamma)
Convert_316:=FourPi/3*sqrt(2/35);   # multiplies Y316 to get [QxQxQ]_(L=6)

# The following specify, in internal format, expressions for the
# quadrupole operator, in two-SU(1,1) basis and single-SU(1,1) basis resp.

quad_op:=[ [Convert_112, [Radial_b,SpHarm_112]] ]:
#quad_op_sq:=[ [Convert_112, [Radial_b_sq,SpHarm_112]] ]:  # not needed??


# The following are functions through which lambda can be determined
# as a function of seniority v.
# (The actual function used is stored in the variable glb_lam_fun,
# and this is only used in the function RepXspace_Twin).
# These functions return lambda_v-lambda_0 which must be an integer.
# And we must have lambda_v-lambda_{v+1} = -1,0 or +1.
# Anything else will cause problems!

# First, completely fixed lambda

lambda_fix_fun:=proc(v::nonnegint)
  0
end;

# Next, the original harmonic oscillator (useful for near-spherical nucleii)

lambda_sho_fun:=proc(v::nonnegint)
  v
end;

# Then the parity values (most useful for well-deformed nucleii).

lambda_acm_fun:=proc(v::nonnegint)
  irem(v,2) 
end;

# A little mixture

lambda_jig_fun:=proc(v::nonnegint)
  if v=0 then
    0
  else
    2-irem(v,2)
  fi
end;


# The following determines functions used for displaying
# the transition rates and amplitudes.  (could avoid using CG)
# (Maple doesn't allow me to specify a delimiting fourth argument $ here!)
# The procedure int_amp_mul, when invoked using ACM_set_transition,
# reproduces what I originally coded in acm13.mpl (described in acm_code04).

quad_rat_fun:=proc(Li,Lf,Mel)
  Mel^2*dimSO3(Lf)/dimSO3(Li)
end;

mel_rat_fun:=proc(Li,Lf,Mel)
  Mel^2*dimSO3(Lf)
end;

unit_rat_fun:=proc(Li,Lf,Mel)
  Mel^2
end;

# We now calculate the Li and Lf dependent scalings
# c.f. eqns. (64) & (65) of [Rowe2004].

quad_amp_fun:=proc(Li,Lf,Mel)
  Mel*gen_amp_mul(Li,Lf,2)
end;

gen_amp_mul:=proc(Li,Lf,Lt,$)
  if Li=Lf then CG(Lf,Lf,Lt,0,Lf,Lf)
  else sqrt(2*Lf+1)
  fi:
end;

int_amp_fun:=proc(Li,Lf,Mel)
  Mel*CG(Li,Li,2,Lf-Li,Lf,Lf)
end;

mel_amp_fun:=proc(Li,Lf,Mel)
  Mel*sqrt(2*Lf+1)
end;

unit_amp_fun:=proc(Li,Lf,Mel)
  Mel
end;

# The following determine how the transition rates and amplitudes
# of the quadrupole operator are displayed by the procedure Show_Rats.

quad_rat_fmat:="  B(E2: %s(%s) -> %s(%s)) = %s":
quad_rat_flist:="  B(E2: %s(#) -> %s(%s)) = %s":
quad_amp_fmat:="    (amplitude = %s)":
quad_amp_flist:="          (Amplitudes = %s)":

# These can be used for other operators, if the user hasn't anything else.

def_rat_fmat:="  Rate( %s(%s) -> %s(%s) ) = %s":
def_rat_flist:="  Rates( %s(#) -> %s(%s) ) = %s":
def_amp_fmat:="    (amplitude = %s)":
def_amp_flist:="          (Amplitudes = %s)":

# And correspondingly for the Show_Mels procedure.

def_mel_fmat:="  ME( %s(%s) -> %s(%s) ) = %s":
def_mel_flist:="  ME( %s(#) -> %s(%s) ) = %s":


# The following defines a table wherein the SO(5) Clebsch-Gordon
# coefficients will be stored in memory. This table is intially empty.
# The table is loaded from external files, as required.
# For a particular (v1,v2,a2,L2,v3), this is done by calling
# load_CG_table(v1,v2,a2,L2,v3).
# When present, the SO(5)>SO(3) reduced CG coefficient is given by
# CG_coeffs[v1,v2,a2,L2,v3][a1,L1,a3,L3].

CG_coeffs:=table():

# To examine which (v1,v2,a2,L2,v3) have been loaded, we can use:
#   indices(CG_coeffs);
# Intially, of course, this table will be empty.


###########################################################################
####---------------------- Global Parameters --------------------------####
###########################################################################

# The SO(5)>SO(3) CG coefficients are initially obtained from external files
# The following determines the directory containing the
# SO(5)>SO(3) CG coefficients.
# This should be specified in the acm_user.mpl file.
# Here we give a sample definition.

SO5CG_directory:="/home/moi/maple/acm/so5_cgs/":

# The following line would test the above directory (it is used by
# the procedure SO5CG_filename), and, somewhat, the data therein
#     readdata( SO5CG_filename(0,3,3),float);  # test
# (it should return a list of four values, each 1.0):


# The data that is produced by the main procedures, is displayed
# according to the values of various global parameters.
# These are listed here, along with some initial values.
# However, it is intended that all of these values will be reset using 
# the user configuable file acm_user.mpl .
# The values here should not be changed directly, but by using the
# routines below.

glb_eig_sft:=1.0:
glb_rat_sft:=1.0:
glb_amp_sft:=1.0:

glb_rel_pre:=2:
glb_rel_wid:=7:
glb_low_pre:=4:

glb_eig_num:=4:
glb_rat_num:=4:

glb_eig_fit:=4.0:
glb_eig_L:=2:
glb_eig_idx:=1:
glb_rat_fit:=10.0:
glb_rat_L1:=2:
glb_rat_L2:=0:
glb_rat_1dx:=1:
glb_rat_2dx:=1:

glb_lam_fun:=lambda_acm_fun:
glb_rat_TRop:=quad_op:
glb_rat_fun:=quad_rat_fun:
glb_amp_fun:=quad_amp_fun:
glb_amp_sft_fun:=sqrt:
glb_rat_format1:=quad_rat_fmat:
glb_rat_format2:=quad_rat_flist:
glb_amp_format1:=quad_amp_fmat:
glb_amp_format2:=quad_amp_flist:
glb_mel_format1:=def_mel_fmat:
glb_mel_format2:=def_mel_flist:

glb_rat_lst:=[]:
glb_amp_flg:=4:

# The ACM_set_scales() procedure may be used to alter the first 3 values.
# The ACM_set_output() procedure may be used to alter the next 3 values.
# The ACM_set_listln() procedure may be used to alter the next 2 values.
# The ACM_set_eig_fit() and ACM_set_rat_fit() procedures may be used to
# alter the next 8 values.
# ACM_set_rat_lst() and ACM_add_rat_lst() affect glb_rat_lst.
# ACM_set_amp_show() affects the last value.


# ACM_set_scales() sets the values of glb_eig_sft, glb_rat_sft, glb_amp_sft
# which are used to divide the values of the (relative) eigenenergies
# and transition rates in the procedure ACM_Scale()
# via the procedures Show_Eigs() and Show_Rats().
# Note that the scaling factor glb_amp_sft is obtained from glb_rat_sft
# using the procedure given by glb_amp_sft_fun.
# This procedure can be changed using ACM_set_amp_fun.
# ACM_set_scales() displays the values of glb_eig_sft, glb_rat_sft, glb_amp_sft.

ACM_set_scales:=proc(eig_sft::numeric,rat_sft::numeric,
                      show::integer:=1,$)
      global glb_eig_sft,glb_rat_sft,glb_amp_sft,glb_amp_sft_fun;

  if nargs>0 then
    glb_eig_sft:=evalf(eig_sft);
  fi:
  if nargs>1 then
    glb_rat_sft:=evalf(rat_sft);
  fi:

  glb_amp_sft:=glb_amp_sft_fun(glb_rat_sft);  #default is square root of above

  if show>0 then
    printf("Relative energies to be multiplied by %f.\n",evalf(1/glb_eig_sft));
    printf("Transition rates to be multiplied by %f.\n",evalf(1/glb_rat_sft));
    printf("Amplitudes to be multiplied by %f.\n\n",evalf(1/glb_amp_sft));
  fi:

  return NULL;
end;

ACM_show_scales:=proc(show::integer:=1,$)
      global glb_eig_sft,glb_rat_sft,glb_amp_sft;

  if show>0 then
    printf("Relative energies to be multiplied by %f.\n",evalf(1/glb_eig_sft));
    printf("Transition rates to be multiplied by %f.\n",evalf(1/glb_rat_sft));
    printf("Amplitudes to be multiplied by %f.\n\n",evalf(1/glb_amp_sft));
  fi:

  return [glb_eig_sft,glb_rat_sft,glb_amp_sft]:
end;

ACM_set_sft_fun:=proc(amp_fun::procedure:=glb_amp_sft_fun,
                         show::integer:=1,$)
    global glb_amp_sft_fun;

  glb_amp_sft_fun:=amp_fun;

  if show>0 then
      printf("Transition amplitudes scaling factor calculated"
                 " using the procedure: \"%a\",\n",
                    glb_amp_sft_fun):
  fi:
end;

# ACM_set_output() sets the values of
# glb_rel_wid, glb_low_pre which are used to format the values output
# by the routines ACM_Scale() and ACM_Adapt(), via
# Show_Eigs() and Show_Rats().

ACM_set_output:=proc(rel_pre::nonnegint,rel_wid::nonnegint,low_pre::nonnegint,
                      show::integer:=1,$)
    global glb_low_pre,glb_rel_wid,glb_rel_pre,glb_eig_num,glb_rat_num;

  if nargs>0 then
    glb_rel_pre:=rel_pre;
  fi:
  if nargs>1 then
    glb_rel_wid:=rel_wid;
  fi:
  if nargs>2 then
    glb_low_pre:=low_pre;
  fi:

  if show>0 then
  printf("%d decimal places for each relative eigenvalue,\n",glb_rel_pre);
  printf("%d total digits for each relative eigenvalue,\n",glb_rel_wid);
  printf("%d decimal places for lowest (absolute) eigenvalue.\n",glb_low_pre);
  fi:

  return NULL;
end;

# ACM_set_listln() sets the values of glb_eig_num, glb_rat_num,
# which determine lengths of displayed lists within
# the routines ACM_Scale() and ACM_Adapt(), via
# Show_Eigs() and Show_Rats().

ACM_set_listln:=proc(eig_num::nonnegint,rat_num::nonnegint,
                      show::integer:=1,$)
      global glb_eig_num,glb_rat_num;

  if nargs>0 then
    glb_eig_num:=eig_num;
  fi:
  if nargs>1 then
    glb_rat_num:=rat_num;
  fi:

  if show>0 then
  printf("Display lowest %d eigenvalue(s) at each L,\n",glb_eig_num);
  printf("Display lowest %d transition rate(s) in each list.\n",glb_rat_num);
  fi:

  return NULL;
end;

# ACM_set_eig_fit() sets the values of
#      glb_eig_fit, glb_eig_L, glb_eig_idx,
# which are used to determine a factor which are used to scale all the
# energies for output. This factor is determined such that the energy
# of the (glb_rat_1dx)th state of AM glb_rat_L1 comes out to be glb_eig_fit.

ACM_set_eig_fit:=proc(eig_fit::numeric, eig_L::nonnegint, eig_idx::posint,
                                                  show::integer:=1,$)
      global glb_eig_fit, glb_eig_L, glb_eig_idx;

  if nargs>0 then
    glb_eig_fit:=eig_fit;
  fi:
  if nargs>1 then
    glb_eig_L:=eig_L;
  fi:
  if nargs>2 then
    glb_eig_idx:=eig_idx;
  fi:

  if show>0 then
    printf("Relative eigenvalue of %d(%d) state set to %f\n",
                                       glb_eig_L,glb_eig_idx,glb_eig_fit);
  fi:

  return NULL;
end;

# In a similar way, ACM_set_rat_fit() sets the values of
#      glb_rat_fit, glb_rat_L1, glb_rat_1dx, glb_rat_L2, glb_rat_2dx:
# which are used to scale the transition rates output.

ACM_set_rat_fit:=proc(rat_fit::numeric, rat_L1::nonnegint, rat_L2::nonnegint,
                       rat_1dx::posint, rat_2dx::posint, show::integer:=1,$)
      global glb_rat_fit, glb_rat_L1, glb_rat_1dx, glb_rat_L2, glb_rat_2dx:

  if nargs>0 then
    glb_rat_fit:=rat_fit;
  fi:
  if nargs>1 then
    glb_rat_L1:=rat_L1;
  fi:
  if nargs>2 then
    glb_rat_L2:=rat_L2;
  fi:
  if nargs>3 then
    glb_rat_1dx:=rat_1dx;
  else
    glb_rat_1dx:=1;
  fi:
  if nargs>4 then
    glb_rat_2dx:=rat_2dx;
  else
    glb_rat_2dx:=1;
  fi:

  if show>0 then
    printf("Transition rate B(E2: %d(%d) -> %d(%d)) set to %f\n",
           glb_rat_L1,glb_rat_1dx,glb_rat_L2,glb_rat_2dx,glb_rat_fit);
  fi:

  return NULL;
end;


# The following three functions change or display the lists of transition
# rates and amplitudes that are flagged for display.
# The functions ACM_set_rat_lst(), ACM_add_rat_lst(), and ACM_show_rat_lst()
# respectively set, append to, and display the current list.

ACM_set_rat_lst:=proc(rat_lst::list(list(integer)),$)
    global glb_rat_lst;

  glb_rat_lst:=[];
  ACM_add_rat_lst(rat_lst):
end;

ACM_add_rat_lst:=proc(rat_lst::list(list(integer)),$)
    local rat_ent;
    global glb_rat_lst;

  for rat_ent in rat_lst do
    if nops(rat_ent)<2 or nops(rat_ent)>5 then
      printf("  Bad transition rate specification: %a\n",rat_ent):
    else
      glb_rat_lst:=[op(glb_rat_lst),rat_ent]:
    fi:
  od:

  return nops(glb_rat_lst);
end;

#

ACM_show_rat_lst:=proc(show::integer:=1,$)
      local rate_ent,rat_format4,rat_format5;
      global glb_rat_lst,glb_rat_format1,glb_rat_format2;

  if show>0 then

    rat_format4:=sprintf(glb_rat_format1,"%d","%d","%d","%d","*"):

    printf("Following transition rates set to be displayed:\n"):
    for rate_ent in glb_rat_lst do

      if nops(rate_ent)=4 or (nops(rate_ent)=5 and rate_ent[5]=0) then
        printf(rat_format4,
             rate_ent[1], rate_ent[3], rate_ent[2], rate_ent[4]):
          printf("\n"):
      elif nops(rate_ent)=5 then
          rat_format5:=sprintf(glb_rat_format1,
                          "%d%+dk","%d","%d%+dk","%d","*"):
          printf(rat_format5, rate_ent[1], rate_ent[5], rate_ent[3],
                              rate_ent[2], rate_ent[5], rate_ent[4]):
        printf("\n"):
      elif nops(rate_ent)=3 then
          rat_format5:=sprintf(glb_rat_format1,
                          "%d","j_i","%d","%d","*"):
          printf(rat_format5, rate_ent[1], rate_ent[2], rate_ent[3]):
        printf("\n"):
      elif nops(rate_ent)=2 then
          rat_format5:=sprintf(glb_rat_format1,
                          "%d","j_i","%d","j_f","*"):
          printf(rat_format5, rate_ent[1], rate_ent[2]):
        printf("\n"):
      fi

    od
  fi:

  return glb_rat_lst;
end;


# The procedure ACM_amp_show(n) indicates for which elements of
# the list of transition rates, are the amplitudes to be displayed
# as well. 6 for none, otherwise only those designators with at
# least n parameters specified.

ACM_set_amp_show:=proc(n::nonnegint,show::integer:=1,$)
  global glb_amp_flg;
  glb_amp_flg:=n:
  if show>0 then
    if n>5 then
      printf("Transition amplitudes not to be displayed.\n"):
    elif n<3 then
      printf("Transition amplitudes displayed for all designators.\n"):
    else
      printf("Transition amplitudes displayed"
                       " only for %d+ element designators.\n",n):
    fi
  fi:

  return NULL;
end;


# The following specifies the operator with respect to which
# transition rates (and amplitudes) are calculated in the
# procedures ACM_Scale and ACM_Adapt.
# (There seems to be nothing sensible to test the 2nd and 3rd
#   argument types against.)

ACM_set_transition:=proc(TR_op::list(list):=glb_rat_TRop,
                         TR_amp_fun::procedure:=glb_amp_fun,
                         TR_rat_fun::procedure:=glb_rat_fun,
                         show::integer:=1,$)
    global glb_rat_TRop,glb_amp_fun,glb_rat_fun;

  glb_rat_TRop:=TR_op;
  glb_amp_fun:=TR_amp_fun;
  glb_rat_fun:=TR_rat_fun;

  if show>0 then
      printf("In ACM_Scale and ACM_Adapt, transition rates "
                 "now calculated for the operator:\n  %a,\n",
                    glb_rat_TRop):
      printf("the transition amplitudes are calculated"
                 " using the procedure: \"%a\",\n",
                    glb_amp_fun):
      printf("and the transition rates are calculated"
                 " using the procedure: \"%a\".\n",
                    glb_rat_fun):
  fi:
end;


# The following specifies the format used to output transition rates
# and amplitudes in the procedure Show_Rats (which is called by
# ACM_Scale and ACM_Adapt. Also set here is op_L which specifies
# the AM of the tensor operator - this is required to correctly
# calculate the transition amplitudes (but does not affect the
# transition rates).

ACM_set_rat_format:=proc(rat_format1::string,rat_format2::string,
                         amp_format1::string,amp_format2::string,
                         show::integer:=1,$)
      global glb_rat_format1,glb_rat_format2,
             glb_amp_format1,glb_amp_format2:

  if nargs>0 then
    glb_rat_format1:=rat_format1;
  fi:
  if nargs>1 then
    glb_rat_format2:=rat_format2;
  fi:
  if nargs>2 then
    glb_amp_format1:=amp_format1;
  fi:
  if nargs>3 then
    glb_amp_format2:=amp_format2;
  fi:

  if show>0 then
    printf("Single transition rates output using formats:\n  ");
    printf(glb_rat_format1,"Li","ji","Lf","jf","*");
    printf("\nTransition rate lists output using format (ji varies):\n  ");
    printf(glb_rat_format2,"Li","Lf","jf","[*]");

    printf("\nSingle amplitudes output using format:\n  ");
    printf(glb_amp_format1,"*");
    printf("\nAmplitude lists output using format:\n  ");
    printf(glb_amp_format2,"[*]");
  fi:

  return NULL;
end;


# Specifies which basis to use, i.e., how lambda varies with v.

ACM_set_basis_type:=proc(choice::nonnegint, abeta0::numeric:=0.0,
                                                  show::integer:=1,$)
  local new_fun:
  global glb_lam_fun, lambda_fix_fun, lambda_sho_fun,
                          lambda_acm_fun, lambda_jig_fun:
  if choice=0 then
    ACM_set_lambda_fun(lambda_fix_fun,0):
    if show>0 then
      printf("Using the constant lambda basis.\n"):
    fi
  elif choice=1 then
    ACM_set_lambda_fun(lambda_sho_fun,0):
    if show>0 then
      printf("Using the harmonic oscillator basis.\n"):
    fi
  elif choice=2 then
    ACM_set_lambda_fun(lambda_acm_fun,0):
    if show>0 then
      printf("Using the ACM parity basis.\n"):
    fi
  elif choice=3 then
    new_fun:=lambda_davi_fun(abeta0^4):
    ACM_set_lambda_fun(new_fun,0):
    if show>0 then
      printf("Using integer Davidson basis for potential with "
             "minimum at %a.\n",abeta0):
    fi
  else
    error "There is no basis %1 defined!", choice:
  fi:

  return NULL;
end;

ACM_set_lambda_fun:=proc(lambda_fun::procedure, show::integer:=1,$)
    global glb_lam_fun;
  glb_lam_fun:=lambda_fun:

  if show>0 then
      printf("lambda values calculated from v"
                 " using the procedure: \"%a\",\n",
                    glb_lam_fun):
  fi:
end;

ACM_show_lambda_fun:=proc(vmin::nonnegint:=0,vmax::nonnegint:=10)
  global glb_lam_fun;
  [seq(glb_lam_fun(v),v=vmin..vmax)]:
end;


###########################################################################
####---------- Procedures for Dimensions and State labels -------------####
###########################################################################

# The following two give the dimensions of SO(3) and SO(5)
# irreducible representations (symmetric).
# The third gives the total number of SO(3) irreps
# (some possibly equivalent) in the SO(5) irrep (v).
# The fourth sums these over the range v_min...v_max.

dimSO3:=(L) -> 2*L+1:                      # SO(3) irrep dimension
dimSO5:=(v) -> (v+1)*(v+2)*(2*v+3)/6:      # SO(5) irrep dimension
dimSO5r3_allL:=(v)                      # irrep count: fixed v, all L
    -> iquo( v*(v+3), 6 ) + 1:
dimSO5r3_rngVallL:=(v_min,v_max)        # irrep count: range v, all L
    -> add(dimSO5r3_allL(v),v=v_min..v_max):
 

# We also need similar functions when L is fixed, or takes a range.
# These make use of the function dimSO5r3(v,L) which gives the
# multiplicity of the states with a given seniority v and angular
# momentum L. This is then the maximum value of the "missing" label
# alpha (the minimum value is 1).

# This formula is given by CRW, eqn. (A1).

dimSO5r3:=proc(v::integer,L::integer)
  local b,d;

  if v<0 or L<0 or L>2*v then
    0:
  else
    b:=(L+3*irem(L,2))/2;
    if v>=b then
      d:=1+iquo(v-b,3);
    else
      d:=0;
    fi:
    if v>=L-2 then
      d:=d-iquo(v-L+2,3);
    fi:
    d:
  fi:
end:

# TAW obtained the following alternative expression.

#alpha_max:=proc(v::integer,L::integer)
#  local q,p,parity;
#  if v<0 or L<0 or L>2*v then
#    0;
#  else
#    parity:=irem(L,2);
#    p:=min(L-3*parity,2*v-L-3*parity) / 2:
#    q:=irem(v+L,3):
#    if p<q then
#      0:
#    else
#      1+iquo(p-q,3):
#    fi:
#  fi:
#end:


dimSO5r3_rngL:=(v,L_min,L_max)                # irrep count: fix v, range L
    -> add( dimSO5r3(v,j), j=L_min..L_max):
dimSO5r3_rngV:=(v_min,v_max,L)            # irrep count: range v, fix L
    -> add( dimSO5r3(i,L),i=v_min..v_max):
dimSO5r3_rngVrngL:=(v_min,v_max,L_min,L_max)  # irrep count: range v and L
    -> add(dimSO5r3_rngL(i,L_min,L_max),i=v_min..v_max):
dimSO5r3_rngVvarL:=(v_min,v_max,L,L_max)
    -> `if`(nargs>3,dimSO5r3_rngVrngL(args),
                    dimSO5r3_rngV(args)):

# We now specify functions which list the labels that correspond
# to the previous dimension formulae.
# These labels are triples [v,alpha,L], where v labels the senioirty
# of an SO(5) irrep, and L labels the SO(3) irreps contained within.
# alpha is a "missing label" which serves to distinguish different
# SO(3) irreps in an SO(5) irrep of seniority v.
# Note that the quantum label M is not included in this list of states
# (it would vary over 2L+1 values: usually -L,-L+1,-L+2,...,L).

# The names of the following functions correspond to those above with
# "dimSO5r3" replaced by "lbsSO5r3".
# The function lbsSO5r3_allL(v) produces a list of pairs [alpha,L]
# that label the SO3 irreps in the SO5 irrep of seniority v.
# All the other functions produce a list of triples [v,alpha,L]
# with v and L within the ranges specified.
# L varies slowest, then v, with the missing label alpha varying fastest.

# When setting up bases for the ACM, later functions use
# lbsSO5r3_rngVvarL(v0,v1,L0,L1), which gives all states in ranges
# for both v and L. The last argument L1 may be omitted whence only
# a single values of L=L0 is used.
# With the L label varying slowest, matrices built using such single L
# cases are readily joined together to give those resulting from a
# range L in L0..L1.

lbsSO5r3_allL:=proc(v::integer)
    [seq(seq([a,LL],a=1..dimSO5r3(v,LL)),LL=0..2*v)]:
end:

lbsSO5r3_rngVallL:=proc(v_min::integer,v_max::integer)
  if v_min<0 or v_min>v_max then
    error("Seniority range invalid");
  else
    [seq(seq(seq([u,a,LL],a=1..dimSO5r3(u,LL)),u=v_min..v_max),
                                                      LL=0..2*v_max)]:
  fi:
end:

# following is used by routines to calculate QixQxQ etc.

lbsSO5r3_rngL:=proc(v::integer,L_min::integer,L_max::integer)
  if v<0 or L_min>L_max then
    error("Parameter range invalid");
  else
    [seq(seq([v,a,LL],a=1..dimSO5r3(v,LL)),LL=L_min..L_max)]:
  fi:
end:

lbsSO5r3_rngV:=proc(v_min::integer,v_max::integer,L::integer)
  if v_min<0 or v_min>v_max then
    error("Seniority range invalid");
  else
    [seq(seq([u,a,L],a=1..dimSO5r3(u,L)),u=v_min..v_max)]:
  fi:
end:

lbsSO5r3_rngVrngL:=proc(v_min::integer,v_max::integer,
                            L_min::integer,L_max::integer)
  if v_min<0 or v_min>v_max or L_min>L_max then
    error("Seniority range invalid");
  else
    [seq(seq(seq([u,a,LL],a=1..dimSO5r3(u,LL)),u=v_min..v_max),
                                  LL=L_min..L_max)]:
  fi:
end:

# The following supersedes the above two - the final (4th) argument may
# be omitted, and if so, only the single value of L is used.

lbsSO5r3_rngVvarL:=proc(v_min::integer,v_max::integer,
                                       L::integer,L_max::integer,$)
  if v_min<0 or v_min>v_max then
    error("Seniority range invalid");
  elif nargs>3 then
    [seq(seq(seq([u,a,LL],a=1..dimSO5r3(u,LL)),u=v_min..v_max),
                                  LL=L..L_max)]:
  else
    [seq(seq([u,a,L],a=1..dimSO5r3(u,L)),u=v_min..v_max)]:
  fi:
end:


###########################################################################
#### SO(5) Clebsch-Gordon coefficients and reps of Spherical Operators ####
###########################################################################

# The following procedure determines where to locate the individual files
# that contain the SO(5)>SO(3) CG coefficients (see above for a test).
# The names of these files are of the form SO5CG_v1_v2-a2-L2_v3.
# The are assumed to lie in directories named SO5CG_v1_v2_v3 which
# themselves lie in directories
# named "v2=1/", "v2=2/", "v2=3/", etc.,
# each a subdirectory of the directory specified in the
# variable SO5CG_directory (see above).

#SO5CG_filename:=proc(v1::integer,v2::integer,v3::integer)
#    cat(SO5CG_directory,"v2=",v2,"/SO5CG_",v1,"_",v2,"_",v3);
#end:
 
SO5CG_filename:=proc(v1::nonnegint,
                     v2::nonnegint,a2::posint,L2::nonnegint,
                     v3::nonnegint)
    cat(SO5CG_directory,"v2=",v2,"/SO5CG_",v1,"_",v2,"_",v3,
    "/SO5CG_",v1,"_",v2,"-",a2,"-",L2,"_",v3);
end:

 
# For v1,v2,a2,L2,v3 (this quintet labels a CG table),
# CG_labels(v1,v2,a2,L2,v3) lists all quartets [alpha1,L1,alpha3,L3],
# for which |L1-L2| <= L3 <= L1+L2.
# (This list is thus smaller than the direct product of the sets obtained
# from two calls to sph_labels).
# The ordering (and length) of the list accords with the order of the
# CG coefficients in the files produced by TAW
# (these files are called SO5CG_v1_v2-a2-L2_v3 in his directories,
# and listed on his webpages (these pages are yet to be updated!!) ):
# therein the "3" labels vary fastest, the "1" labels next fastest.
# Thus we use the output of this routine to correctly label the data
# from the data file.


CG_labels:=proc(v1::nonnegint,
                v2::nonnegint,a2::posint,L2::nonnegint,
                v3::nonnegint)
  local L1,L3,a1,a3,label_list;
  label_list:=[]:

  for L1 from 0 to 2*v1 do
  for a1 to dimSO5r3(v1,L1) do
    for L3 from abs(L1-L2) to min(L1+L2,2*v3) do
    for a3 to dimSO5r3(v3,L3) do
        label_list:=[op(label_list),[a1,L1,a3,L3]]:
    od: od:
  od: od:

  label_list:
end:


# The following two functions, which read the SO(5) CG coefficients
# for a particular (v1,v2,a2,L2,v3) from the relevant datafile,
# are almost the same (the first returns the data,
# the second prints it out). They are not used subsequently.
 
get_CG_file:=proc(v1::nonnegint,
                  v2::nonnegint,a2::posint,L2::nonnegint,
                  v3::nonnegint)
  local CG_data,CG_list:

  if v1>v2+v3 or v2>v3+v1 or v3>v1+v2 or type(v1+v2+v3,odd)
      or v3<v1        # data is obtained from v3>v1 cases
      or a2>dimSO5r3(v2,L2) then
    error "No CG file for these parameters!":
  fi:

  CG_data:=readdata( SO5CG_filename(v1,v2,a2,L2,v3), float):
  CG_list:=CG_labels(v1,v2,a2,L2,v3):
  [CG_data,CG_list]:
end:

show_CG_file:=proc(v1::nonnegint,
                   v2::nonnegint,a2::posint,L2::nonnegint,
                   v3::nonnegint)
  local cg_table,count,i:

  cg_table:=get_CG_file(v1,v2,a2,L2,v3):
  count:=nops(cg_table[1]):
  if count=1 then
    print(`This file contains ` || count || ` CG coefficient`):
  else
    print(`This file contains ` || count || ` CG coefficients`):
  fi:

  for i to count do
    print(cg_table[2][i],cg_table[1][i]):
  od:
end:


# The following function loads all the SO(5) CG coefficients for a
# particular (v1,v2,a2,L2,v3) from the datafile as determined above.
# Subsequent attemps to load the same values will be ignored with
# no message.
# The data is loaded into the CG_coeffs table, which was initialised above.
# The ordering of the coefficients in the file is important, and assumed
# to be the same as that listed by the routine CG_labels() above.
# Only the first entry (the coefficient) from each line of the file
# is read (we could thus make the data files smaller by only giving
# that one value in each line).

# Note that if v1>v3, then the data is instead loaded for (v3,v2,a2,L2,v1)
# because the CG coefficients for (v1,v2,a2,L2,v3) are easily obtained
# from the former.

# Also note that no checking is done here on the correct ranges
# of the arguments (this is left to the functions that call this).
# We should restrict to
#   v1<=v2+v3 and v2<=v3+v1 and v3<=v1+v2 and type(v1+v2+v3,odd)
#      and a2<=dimSO5r3(v2,L2).

load_CG_table:=proc(v1::nonnegint,
                    v2::nonnegint,a2::posint,L2::nonnegint,
                    v3::nonnegint)
  local CG_data,CG_list,vt1,vt3;
  global CG_coeffs;

  if v1>v3 then
    vt1:=v3: vt3:=v1:
  else
    vt1:=v1: vt3:=v3:
  fi:
  if evalb([vt1,v2,a2,L2,vt3] in [indices(CG_coeffs)] ) then
    return:
  fi:

  CG_data:=readdata( SO5CG_filename(vt1,v2,a2,L2,vt3), float):
  CG_list:=CG_labels(vt1,v2,a2,L2,vt3):

  CG_coeffs[vt1,v2,a2,L2,vt3]:=table([seq( (op(CG_list[i]))=CG_data[i],
                                 i=1..nops(CG_list) )]);
end:


# The following returns the SO(5) CG coefficient
# (v1,a1,L1;v2,a2,L2||v3,a3,L3)    [no renormalisation required],
# loading data from the (v1,v2,a2,L2,v3) datafile if not already loaded.
# Note that if v1>v3, then the data from the (v3,v2,a2,L2,v1) dataset is used,
# and a factor is included (see (4.164) of [RowanWood]).
# (Possibly write a faster version which does no testing of indices).

CG_SO5r3:=proc(v1::nonnegint,a1::posint,L1::nonnegint,
               v2::nonnegint,a2::posint,L2::nonnegint,
               v3::nonnegint,a3::posint,L3::nonnegint)
  global CG_coeffs;
  if v1+v2<v3 or v1+v3<v2 or v2+v3<v1 or
     L1+L2<L3 or L1+L3<L2 or L2+L3<L1 or
     a1>dimSO5r3(v1,L1) or a2>dimSO5r3(v2,L2)
                        or a3>dimSO5r3(v3,L3) or type(v1+v2+v3,odd) then
        0;
  else
     load_CG_table(v1,v2,a2,L2,v3);
     if v1<=v3 then
       CG_coeffs[v1,v2,a2,L2,v3][a1,L1,a3,L3]:
     else
       CG_coeffs[v3,v2,a2,L2,v1][a3,L3,a1,L1]
         * (-1)^(L3+L2-L1)
         * sqrt( dimSO5(v3) * dimSO3(L1) / dimSO5(v1) / dimSO3(L3) ):
     fi:
  fi:
end:


# The following gives the doubly reduced matrix element <u|||w|||v>*4*Pi
# for SO(5) spherical harmonic tensor operators.
# (think of as v==v_i, u==v_f). See Eqn. (4.175) of [RowanWood].
# The return value is algebraic and exact (and probably a surd).
# The normalisation is now correct (an extra factor of sqrt(2)
# was added 17/3/2010!).

ME_SO5red:=proc(u::nonnegint,w::nonnegint,v::nonnegint)
  local sigma,halfsigma;
  if u+v<w or u+w<v or v+w<u or type(u+v+w,odd) then
    RETURN(0);
  fi:

  sigma:=(v+w+u); halfsigma:=sigma/2;
  (halfsigma+1)! / (halfsigma-u)! / (halfsigma-v)! / (halfsigma-w)!
    * sqrt( (2*v+3) * (2*w+3) * (sigma+4) / (u+2) / (u+1)
            * (sigma-2*u+1)! * (sigma-2*w+1)! * (sigma-2*v+1)! / (sigma+3)! );
end:

# The following nine functions are instances of the above, with
# different normalisations: they provide SO(5) (doubly) reduced
# matrix elements for Q and [QxQ]_(v=2) and [QxQxQ]_(v=3)
# From ME_SO5red, they may be obtained using
#   Q=4*Pi/sqrt(15) * Y_{112} and
#   [QxQ]_(24)=4*Pi*sqrt(2/105) * Y_{214} and
#   [QxQxQ]_(36)=4*Pi*sqrt(2/315) * Y_{310},
# (we need to take care with different signs for
#       [QxQ]_(22)=-4*Pi*sqrt(2/105) * Y_{212} and
#       [QxQxQ]_(30)=-4*Pi*sqrt(2/315) * Y_{316}).

# These give the expressions in large parentheses in eqs (A2-4),
# and thus after multiplication by a SO(5)>SO(3) CG coefficient,
# and possibly a sign, give SO(3) reduced matrix elements of
# [Q], [QxQ]_{v=2,L=2,4} and [QxQxQ]_{v=3,L=0,6}.

# Later, we use these in procedures to calculate [Q^+ x Q x Q]_{v=3}
# and [Q^- x Q x Q]_{v=3}.

Qred_p1:=(v)->sqrt((v+1)/(2*v+5));
Qred_m1:=(v)->sqrt((v+2)/(2*v+1));

QxQred_p2:=(v)->sqrt((v+1)*(v+2)/(2*v+5)/(2*v+7));
QxQred_0:=(v)->sqrt(6*v*(v+3)/5/(2*v+1)/(2*v+5));
QxQred_m2:=(v)->sqrt((v+1)*(v+2)/(2*v+1)/(2*v-1));

QxQxQred_p3:=(v)->sqrt((v+1)*(v+2)*(v+3)/(2*v+5)/(2*v+7)/(2*v+9));
QxQxQred_p1:=(v)->3*sqrt(v*(v+1)*(v+4)/7/(2*v+1)/(2*v+5)/(2*v+7));
QxQxQred_m1:=(v)->3*sqrt((v-1)*(v+2)*(v+3)/7/(2*v-1)/(2*v+1)/(2*v+5));
QxQxQred_m3:=(v)->sqrt(v*(v+1)*(v+2)/(2*v-3)/(2*v-1)/(2*v+1));



# The following obtains the (adjusted SO(3) reduced) matrix element
#
#      4*Pi
#   ------------- * <v_f,al_f,L_f || Y_{v,al,L} || v_i,al_i,L_i>
#   sqrt (2*L_f+1)
#
# for the SO(5) spherical harmonic Y_{v,al,L}.
# (See eqn. (51) of [RWC2009], and compare with eqn. (49).)
# Note that the case where v_i>v_f is dealt with correctly
# (the jiggery in CG_SO5r3() together with the form of the expression
# in ME_SO5red() correctly reproduces [RowanWood]'s (3.155) ).

# To get things like cos(3g), we use the relation Y310 = (3/4/Pi) cos(3g).
# So to get the matrix element of cos(3g) from this, we must mult by (1/3).

# The return value may be (partly) algebraic, so it may be necessary to
# apply evalf to it.

ME_SO5r3:=proc(v_f::integer,al_f::integer,L_f::integer,
               v::integer,al::integer,L::integer,
               v_i::integer,al_i::integer,L_i::integer)

   CG_SO5r3(v_i,al_i,L_i,v,al,L,v_f,al_f,L_f) * ME_SO5red(v_f,v,v_i):
end;


# Now produce a full matrix of the (adjusted SO(3) reduced) matrix elements
# defined above
# (i.e., the result is still lacking sqrt(2L_f+1)/4/Pi factors).
# The seniorities range from v_min to v_max.
# BEWARE with the AM arguments: the L_max argument is optional
#  - if omitted, the single value of AM is used; and if L_max is given,
# all values between Lmin and Lmax (inclusive) are used.

# The remember option is used, but these stored values are cleared
# each time the (much later) RepXspace_OpLC() routine is invoked.

# The evalf within this procedure ensures that the returned matrix
# has float entries.

# If [v,al,L] doesn't label a spherical harmonic, then a matrix of 0s is
# returned. (Perhaps this latter behaviour should change.)

RepSO5_Y:=proc(v::integer,al::integer,L::integer,
                    v_min::integer,v_max::integer,
                    L_min::integer,L_max::integer,$)
  option remember;
  local states:

  states:=lbsSO5r3_rngVvarL(args[4..-1]): 
  Matrix( nops(states), (i,j)->evalf(
                 ME_SO5r3(op(states[i]),v,al,L,op(states[j])) )):
end:


# As above, but without evalf so that the ME's come out (partially) algebraic.

RepSO5_Y_alg:=proc(v::integer,al::integer,L::integer,
                    v_min::integer,v_max::integer,
                    L_min::integer,L_max::integer,$)
  local states:

  states:=lbsSO5r3_rngVvarL(args[4..-1]): 
  Matrix( nops(states), (i,j)->
                 ME_SO5r3(op(states[i]),v,al,L,op(states[j])) ):
end:



# Forms the diagonal matrix whose entries are (-1)^L*sqrt(2L+1).

RepSO5_sqLdim:=proc(v_min::integer,v_max::integer,
                       L::integer,L_max::integer,$)
  local states:

  # Obtain labels for space.

  states:=lbsSO5r3_rngVvarL(args): 
 
  # now form the diagonal matrix, of sqrt(2L+1)

  Matrix(map(x->evalf(eval((-1)^(x[3])*sqrt(dimSO3(x[3])))),states),
                               shape=diagonal,scan=diagonal);
end:

# Forms the diagonal matrix whose entries are (-1)^L/sqrt(2L+1).

RepSO5_sqLdiv:=proc(v_min::integer,v_max::integer,
                       L::integer,L_max::integer,$)
  local states:

  # Obtain labels for space.

  states:=lbsSO5r3_rngVvarL(args): 
 
  # now form the diagonal matrix, of 1/sqrt(2L+1)

  Matrix(map(x->evalf(eval((-1)^(x[3])/sqrt(dimSO3(x[3])))),states),
                               shape=diagonal,scan=diagonal);
end:


# The next procedure calculates the representation matrix of a product
# of Spherical Harmonics (only certain products are meaningful).
# This product is encoded in ys_op as an array [a1,a2,a3,a4,.....],
# with each element one of the symbolic names
#     SpHarm_112, SpHarm_212, SpHarm_214,
#     SpHarm_310, SpHarm_313, SpHarm_314, SpHarm_316,
#     ...... SpHarm_61C.
# or is a triple [v,alpha,L] which designates the Spherical harmonic
# explicitly, or is one of the diagonal operators
#     SpDiag_sqLdim, SpDiag_sqLdiv.
# Because of the use of RepSO5_Y above, the matrix returned has
# entries of type float.
# The states are those [v',alpha',L'] with v' in the range v_min..v_max.
# BEWARE with the L arguments: the L_max argument is optional - if omitted,
# the single value of L_min is used;
# and if L_max is given, we use all values between L_min and Lmax (inclusive).
# If [v',alpha',L'] does not correspond to an actual spherical harmonic
# (in our table), then a zero matrix is used.
#
# No renormalisation to other operators is performed here (done elsewhere).
# These are for the pure orthonormal spherical harmonics.
# However, to obtain adjusted SO(3)-reduced matrix elements the
# result needs to be multiplied by (4*Pi)^(-n) where n is the
# number of spherical harmonics in the list ys_op
# (n may be obtained using the procedure NumSO5r3_Prod() below,
# and is not necessarily the length of ys_op because this list may
# contain non-harmonics such as SpDiag terms).
# To obtain genuine SO(3)-reduced matrix elements, the result
# needs to be further multiplied by (a matrix of) sqrt(2L_f+1).
#
# Note the use of "copy" in this procedure. This is necessary otherwise
# the remember table for RepSO5_Y gets messed up!

# It might be a good idea also to use remember on this procedure,
# because this might get called repeatedly for a particular Hamiltonian.


RepSO5r3_Prod:=proc(ys_op::list,v_min::integer,v_max::integer,
                                         L_min::integer,L_max::integer,$)
    local i,n,Mats,Mat_product,this_op;
    global SpHarm_Indices,SpHarm_Table;

  n:=nops(ys_op);

  if n=0 then   # require identity matrix
    return Matrix([seq(1,i=1..dimSO5r3_rngVvarL(args[2..-1]))],
                                   scan=diagonal):
  else

    for i from 1 to n do

      if type(ys_op[i],list(nonnegint)) and nops(ys_op[i])=3 then
        this_op:=ys_op[i];
      elif type(ys_op[i],nonnegint) and member(ys_op[i],SpHarm_Indices) then
        this_op:=SpHarm_Table[ys_op[i]];
      elif type(ys_op[i],nonnegint) and ys_op[i]=SpDiag_sqLdim then
        if i=1 then  # cannot make copy because that'd force diag matrix
          Mat_product:=Matrix(RepSO5_sqLdim(args[2..-1]));
        else
          MatrixMatrixMultiply(Mat_product,
                               RepSO5_sqLdim(args[2..-1]),inplace);
        fi:
        next;     # tackle next i in for loop
      elif type(ys_op[i],nonnegint) and ys_op[i]=SpDiag_sqLdiv then
        if i=1 then
          Mat_product:=Matrix(RepSO5_sqLdiv(args[2..-1]));
        else
          MatrixMatrixMultiply(Mat_product,
                               RepSO5_sqLdiv(args[2..-1]),inplace);
        fi:
        next;     # tackle next i in for loop

      else
        error "Invalid SO(5) harmonic designator %1", ys_op[i]:
      fi:


      # Now multiply in the spherical harmonic denoted by this_op.

      if i=1 then
        Mat_product:=copy(RepSO5_Y( op(this_op), args[2..-1]));
      else
        MatrixMatrixMultiply(
                   Mat_product,
                   RepSO5_Y( op(this_op), args[2..-1]),
                   inplace);
      fi:

    od:
  fi:
  
#  To get genuine adjusted reduced matrix elements, we now need to
#  multiply by (4*Pi)^(-n) for n the number of SpHarms in the list

   Mat_product;
end:


# Counts the number of actual spherical harmonics in the list.

NumSO5r3_Prod:=proc(ys_op::list,$)
  local i,ct:

  ct:=0:

  for i from 1 to nops(ys_op) do
    if type(ys_op[i],list(nonnegint)) and nops(ys_op[i])=3 then
        ct:=ct+1:
    elif type(ys_op[i],nonnegint) and member(ys_op[i],SpHarm_Indices) then
        ct:=ct+1:
    fi:
  od:

  ct;
end:



###########################################################################
####--------------- SO(3) Clebsch-Gordon coefficients -----------------####
###########################################################################

# The following function provides the usual SO(3) Clebsch-Gordon
# coefficients CG(J1,m1,J2,m2;J,m). The only place where this might
# get used is where transition amplitudes are calculated from
# transition rates. It's not used elsewhere.
# (Code from DJR)

# Using this, we could, for example get full matrix elements of Y310
# for all states labelled by (v,alpha,L,M) where M is the magnetic
# quantum number, using
#    CG(L_i,M_i,0,0,L_f,M_f)*Y310_ME(v_f,al_f,L_f,v_i,al_i,L_i):

CG:=proc(j1,m1,j2,m2,j3,m3)
  Wigner_3j(j1,j2,j3,m1,m2,-m3)*(-1)^(j1-j2+m3)*sqrt(2*j3+1);
end:

Wigner_3j:=proc(j1,j2,j3,m1,m2,m3)
    local s,ss,f,edmonds;

    if type([args],list(rational)) then
        if abs(m1)>j1 or abs(m2)>j2 or abs(m3)>j3 then RETURN (0) fi;
        if m1+m2+m3 <> 0 then RETURN (0) fi;
        if j3 < abs(j1-j2) or j3 > j1+j2 then RETURN (0) fi;
        ss:=0:
        for s from 0 to j3+m3 do
            if j1-m1-s < 0 or j3+m3-s < 0 or j2-j3+m1+s < 0 then
                 f:=0
           else
               f:=(j1+m1+s)!*(j2+j3-m1-s)!/
                  (s!*(j1-m1-s)!*(j2-j3+m1+s)!*(j3+m3-s)!):
           fi;
           ss:=ss+(-1)**s*f;
       od;
       f:=(j1+j2-j3)!*(j1-m1)!*(j2-m2)!*(j3-m3)!*(j3+m3)!/
          ((j1+j2+j3+1)!*(j1-j2+j3)!*(-j1+j2+j3)!*(j1+m1)!*(j2+m2)!):
       edmonds:=(-1)^(-2*j1-m1-j2-m3)*ss*sqrt(f);
       simplify(edmonds)
   else
       f:=simplify((-1)^(-2*j1-m1-j2-m3)*sqrt((j1+j2-j3)!
          *(j1-m1)!*(j2-m2)!*(j3-m3)!*(j3+m3)!/
          ((j1+j2+j3+1)!*(j1-j2+j3)!*(-j1+j2+j3)!*(j1+m1)!*(j2+m2)!)) *
          add((-1)^s*(j1+m1+s)!*(j2+j3-m1-s)!/
              (s!*(j1-m1-s)!*(j2-j3+m1+s)!*(j3+m3-s)!),s=0..j3+m3));
   fi;
end:


###########################################################################
####----------- Dimensions and labels for Radial and Xspace -----------####
###########################################################################

# The next set of routines deal with representing operators in the radial
# coordinate (beta) space.
# The action is dependent on a parameter lambda (which is sometimes
# regarded as a quantum number).
# Later, we will form a direct product space between these radial states
# and the SO(5) spherical harmonic states. Then, lambda will actually
# be dependent on the seniority v quantum number of the SO(5) labels.


# The next function simply gives the number (dimension) of nu labels
# in the range. (Usually nu_min=0.)
# The one after lists them.

dimRadial:=(nu_min,nu_max) -> nu_max-nu_min+1:

lbsRadial:=proc(nu_min::integer,nu_max::integer)
  if nu_min<0 or nu_min>nu_max then
    error("Radial range invalid");
  else
    [seq(i,i=nu_min..nu_max)];
  fi:
end:


# Give the dimension for the product space nu x sph.
# The last parameter is optional: if not given, L is constant single value.

dimXspace:=proc(nu_min::integer,nu_max::integer,
                   v_min::integer,v_max::integer,L::integer,L_max::integer,$)

    dimRadial(nu_min,nu_max)*dimSO5r3_rngVvarL(args[3..-1]):
end:


# Now list the set of labels for that product space nu x sph.
# Each label is of the form [nu,v,alpha,L].
# The nu label varies fastest, then alpha, then v, and L
# is slowest (as elsewhere).
# The last parameter is optional: if not given, L is constant.

lbsXspace:=proc(nu_min::integer,nu_max::integer,
                   v_min::integer,v_max::integer,L::integer,L_max::integer,$)
    local b_labels,sph_labels;

  b_labels:=lbsRadial(nu_min,nu_max);         # radial labels
  sph_labels:=lbsSO5r3_rngVvarL(args[3..-1]):   # SO(5) labels

  [seq( seq( [nu,op(s)], nu in b_labels), s in sph_labels)];
end:


###########################################################################
####------------------ Representing Radial Operators ------------------####
###########################################################################

# The functions that follow calculate matrix elements
# <R^lambda_mu_f | Op | R^lambda_mu_i> 
# of various operators on the radial space whose elements
# are labelled by mu=0,1,2,3,... .
# For a fixed value of lambda, these elements span the radial space.
# Take lambda>1, else some of these matrix elements are undefined.
# Here, some of the operators link spaces of different lambda.

# Firstly, here are matrix elements of the SU(1,1) operators S0,S+,S-.
# These are as in eqn. (4.94-6) of [RowanWood].

ME_Radial_S0:=proc(lambda::algebraic,mu_f::nonnegint,mu_i::nonnegint)
  if mu_f=mu_i then
    lambda/2 + mu_i;
  else
    0;
  fi;
end:

ME_Radial_Sp:=proc(lambda::algebraic,mu_f::nonnegint,mu_i::nonnegint)
  if mu_f=mu_i+1 then
    sqrt( (lambda + mu_i)*(mu_i+1) );
  else
    0;
  fi;
end:

ME_Radial_Sm:=proc(lambda::algebraic,mu_f::nonnegint,mu_i::nonnegint)
  if mu_f=mu_i-1 then
    sqrt( (lambda + mu_i - 1)*mu_i );
  else
    0;
  fi;
end:


# The following are for products of an even number of beta,
# 1/beta and d/d(beta).
# (Each of these may, in fact, be inferred from the above three.)
# Actually, there is a scale parameter (anorm) knocking around,
# and these procedures actually give the matrix elements of the
# above functions with beta -> (anorm*beta).
# These are as in eqn. (4.97-9) of [RowanWood].

# First, beta^2 (*a^2).

ME_Radial_b2:=proc(lambda::algebraic,mu_f::nonnegint,mu_i::nonnegint)
  if mu_f=mu_i-1 then
    sqrt( (lambda + mu_i - 1)*mu_i );
  elif mu_f=mu_i then
    lambda + 2*mu_i;
  elif mu_f=mu_i+1 then
    sqrt( (lambda + mu_i)*(mu_i+1) );
  else    
    0;
  fi;
end:

# Next, 1/beta^2 (*a^(-2)). It uses the subsequent procedure which
# assumes that mu_f >= mu_i (and restrictions on lambda from Sbm2_ME).

ME_Radial_bm2:=proc(lambda::algebraic,mu_f::nonnegint,mu_i::nonnegint)
  if lambda=-1 then
    error "Singular 1/beta^2 for lambda=1";
  fi:
  if frac(lambda)=0 and (lambda <= -mu_i or lambda <= -mu_f) then
    error "cannot evaluate Gamma function at non-positive integer":
  fi:

  if mu_f>=mu_i then
    ME_Radial_pt(lambda,mu_f,mu_i);
  else
    ME_Radial_pt(lambda,mu_i,mu_f);
  fi:
end:

ME_Radial_pt:=proc(lambda::algebraic,mu_f::nonnegint,mu_i::nonnegint)
  (-1)^(mu_f-mu_i) * sqrt( (factorial(mu_f)*GAMMA(lambda+mu_i))
                               /(factorial(mu_i)*GAMMA(lambda+mu_f)) )
                  / (lambda-1);
end:

# Now it's d^2/d(beta)^2 (*a^(-2)).

ME_Radial_D2b:=proc(lambda::algebraic,mu_f::nonnegint,mu_i::nonnegint)
  local stuff:
  if mu_f=mu_i-1 then
    stuff:=sqrt( (lambda + mu_i - 1)*mu_i );
  elif mu_f=mu_i then
    stuff:=-lambda - 2*mu_i;
  elif mu_f=mu_i+1 then
    stuff:=sqrt( (lambda + mu_i)*(mu_i+1) );
  else    
    stuff:=0;
  fi;

  if mu_f>=mu_i then
    stuff+(lambda-(3/2))*(lambda-(1/2))*ME_Radial_pt(lambda,mu_f,mu_i);
  else
    stuff+(lambda-(3/2))*(lambda-(1/2))*ME_Radial_pt(lambda,mu_i,mu_f);
  fi:
end:

# Next, beta*d/d(beta).  See eqn. (25) of [RWC2009]. It doesn't seem to
# be in [RowanWood], but arises by combining (4.132) & (4.134) therein.

ME_Radial_bDb:=proc(lambda::algebraic,mu_f::nonnegint,mu_i::nonnegint)
  if mu_f=mu_i-1 then
    sqrt( (lambda + mu_i - 1)*mu_i );
  elif mu_f=mu_i then
    -(1/2);
  elif mu_f=mu_i+1 then
    -sqrt( (lambda + mu_i)*(mu_i+1) );
  else    
    0;
  fi;
end:

## The operator 1; changing lambda by -2
#
#ME_Radial_id_ml2:=proc(lambda::algebraic,mu_f::nonnegint,mu_i::nonnegint) 
#  if frac(lambda)=0 and lambda <= -mu_i+2 then
#    error "cannot evaluate Gamma function at non-positive integer":
#  fi:
#
#  if mu_f>mu_i+1 then
#    0
#  elif mu_f=mu_i+1 then
#    sqrt((mu_f)/(lambda+mu_f-2))
#  else
#    (-1)^(mu_f-mu_i)*(lambda-2)
#           *sqrt( (factorial(mu_i)*GAMMA(lambda+mu_f-2))
#                      /(factorial(mu_f)*GAMMA(lambda+mu_i)) )
#  fi:
#end:
#
## The operator 1; changing lambda by +2
#
#ME_Radial_id_pl2:=proc(lambda::algebraic,mu_f::nonnegint,mu_i::nonnegint) 
#  if frac(lambda)=0 and lambda <= -mu_i then
#    error "cannot evaluate Gamma function at non-positive integer":
#  fi:
#
#  if mu_i>mu_f+1 then
#    0
#  elif mu_i=mu_f+1 then
#    sqrt((mu_i)/(lambda+mu_i))
#  else
#    (-1)^(mu_f-mu_i)*lambda
#           *sqrt( (factorial(mu_f)*GAMMA(lambda+mu_i))
#                      /(factorial(mu_i)*GAMMA(lambda+mu_f+2)) )
#  fi:
#end:

# The operator 1; changing lambda by +2*k
# (by swapping mu_f & mu_i and changing lambda -> lambda-2*k,
#  we get the change by -2*k instead, but need a separate function!)

ME_Radial_id_pl:=proc(lambda::algebraic,mu_f::nonnegint,mu_i::nonnegint,
                                  k::nonnegint) 
  if frac(lambda)=0 and lambda <= -mu_i then
    error "cannot evaluate Gamma function at non-positive integer":
  fi:

  if mu_i<=mu_f+k then
    MF_Radial_id_pl(lambda,mu_f,mu_i,k)
            *sqrt( (factorial(mu_f)*GAMMA(lambda+mu_i))
                       /(factorial(mu_i)*GAMMA(lambda+mu_f+2*k)) )
  else
    0
  fi:
end:

# The same operator 1; changing lambda by -2*k

ME_Radial_id_ml:=proc(lambda::algebraic,mu_f::nonnegint,mu_i::nonnegint,
                                  k::nonnegint) 
  if frac(lambda)=0 and lambda <= -mu_f+2*k then
    error "cannot evaluate Gamma function at non-positive integer":
  fi:

  if mu_f<=mu_i+k then
    MF_Radial_id_pl(lambda-2*k,mu_i,mu_f,k)
            *sqrt( (factorial(mu_i)*GAMMA(lambda+mu_f-2*k))
                       /(factorial(mu_f)*GAMMA(lambda+mu_i)) )
  else
    0
  fi:
end:



# The result of the following will be polynomial in lambda, mu and nu.

MF_Radial_id_pl:=proc(lambda::algebraic,mu::nonnegint,nu::nonnegint,
                                  k::nonnegint)
  local res:

  if nu>mu+k then
    return(0):
  fi:

  res:=add( (-1)^j * binomial(k,j) * binomial(k+mu-nu+j-1,k-1)
              * GAMMA(lambda+mu+2*k)/GAMMA(lambda+mu+k+j)
              * GAMMA(mu+j+1)/GAMMA(mu+1),
                                              j=max(0,nu-mu)..k):

  simplify(res)*(-1)^(mu+nu):
end;

# Same result, but done in a different way.

MF_Radial_id_pl2:=proc(lambda::algebraic,mu::nonnegint,nu::nonnegint,
                                   k::nonnegint)
  local res:

  if nu>mu+k then
    return(0):
  fi:

  res:=add( (-1)^j * binomial(k,j) * binomial(2*k+mu-nu-j-1,k+mu-nu)
              * GAMMA(lambda+2*k+mu)/GAMMA(lambda+k-1)
              * GAMMA(lambda+2*k-j-1)/GAMMA(lambda+2*k+mu-j),
                                              j=0..k-1):

  if nu=mu+k then
    res:=res+ (-1)^k * GAMMA(lambda+2*k+mu)/GAMMA(lambda+k+mu):
  fi:

  simplify(res)*(-1)^(mu+nu):
end;




# Now for the matrix elements of operators involving odd
# products of beta and d/d(beta).
# Note that these act between states with different lambda's.
# The initial lambda is that given as the first argument.

# Here's beta; it changes lambda by +1 (eqn. (4.100) of [RowanWood]).

ME_Radial_b_pl:=proc(lambda::algebraic,mu_f::nonnegint,mu_i::nonnegint)
  if mu_f=mu_i-1 then
    sqrt( mu_i );
  elif mu_f=mu_i then
    sqrt(lambda + mu_i);
  else    
    0;
  fi;
end:

# Here's 1/beta; changing lambda by +1 (eqn. (4.105) of [RowanWood]).

ME_Radial_bm_pl:=proc(lambda::algebraic,mu_f::nonnegint,mu_i::nonnegint) 
  if frac(lambda)=0 and lambda <= -mu_i then
    error "cannot evaluate Gamma function at non-positive integer":
  fi:

  if mu_f<mu_i then
    0;
  else
    (-1)^(mu_f-mu_i)*sqrt( (factorial(mu_f)*GAMMA(lambda+mu_i))
                           /(factorial(mu_i)*GAMMA(lambda+mu_f+1)) );
  fi:
end:

# Here's d/d(beta); changing lambda by +1 (eqn. (4.106) of [RowanWood]).

ME_Radial_Db_pl:=proc(lambda::algebraic,mu_f::nonnegint,mu_i::nonnegint) local res:
  if mu_f=mu_i-1 then
    res:=sqrt( mu_i );
  elif mu_f=mu_i then
    res:=-sqrt( lambda + mu_i );
  else    
    res:=0;
  fi;

  if mu_f>=mu_i then
    res:=res+(-1)^(mu_f-mu_i) * (lambda-1/2)
               * sqrt( (factorial(mu_f)*GAMMA(lambda+mu_i))
                          /(factorial(mu_i) * GAMMA(lambda+mu_f+1)) );
  fi:
  res:
end:


# Here's beta; changing lambda by -1 (eqn. (4.103) of [RowanWood]).

ME_Radial_b_ml:=proc(lambda::algebraic,mu_f::nonnegint,mu_i::nonnegint)
  if mu_f=mu_i+1 then
    sqrt( mu_f );
  elif mu_f=mu_i then
    sqrt(lambda + mu_i - 1);
  else    
    0;
  fi;
end:

# Here's 1/beta; changing lambda by -1 (eqn. (4.101) of [RowanWood]).

ME_Radial_bm_ml:=proc(lambda::algebraic,mu_f::nonnegint,mu_i::nonnegint) 
  if frac(lambda)=0 and lambda <= -mu_i then
    error "cannot evaluate Gamma function at non-positive integer":
  fi:

  if mu_f>mu_i then
    0;
  else
    (-1)^(mu_f-mu_i)*sqrt( (factorial(mu_i)*GAMMA(lambda+mu_f-1))
                           /(factorial(mu_f)*GAMMA(lambda+mu_i)) );
  fi:
end:

# Here's d/d(beta); changing lambda by -1 (eqn. (4.102) of [RowanWood]).

ME_Radial_Db_ml:=proc(lambda::algebraic,mu_f::nonnegint,mu_i::nonnegint) local res:
  if mu_f=mu_i+1 then
    res:=-sqrt( mu_f );
  elif mu_f=mu_i then
    res:=sqrt(lambda + mu_i - 1);
  else    
    res:=0;
  fi;

  if mu_f<=mu_i then
    res:=res+(-1)^(mu_f-mu_i) * (3/2-lambda)
               * sqrt( (factorial(mu_i)*GAMMA(lambda+mu_f-1))
                          /(factorial(mu_f)*GAMMA(lambda+mu_i)) );
  fi:
  res:
end:


# The following function calls one of the above to obtain a single
# matrix element in the case of
#   beta^2; 1/beta^2; d^2/d(beta)^2; beta*d/d(beta);
#   beta; 1/beta; d/d(beta)  (each of these three with lambda increasing);
#   beta; 1/beta; d/d(beta)  (each of these three with lambda decreasing);
#   S0; S+; S-.

ME_Radial:=proc(radial_op::integer, anorm::algebraic,
                   lambda::algebraic, lambda_var::integer,
                   mu::nonnegint, nu::nonnegint)
  local MM:

  if radial_op=Radial_b2 and lambda_var=0 then
    ME_Radial_b2(lambda,mu,nu)/anorm^2;
  elif radial_op=Radial_bm2 and lambda_var=0 then
    ME_Radial_bm2(lambda,mu,nu)*anorm^2;
  elif radial_op=Radial_D2b and lambda_var=0 then
    ME_Radial_D2b(lambda,mu,nu)*anorm^2;
  elif radial_op=Radial_bDb and lambda_var=0 then
    ME_Radial_bDb(lambda,mu,nu);        

  elif radial_op=Radial_b and lambda_var=1 then
    ME_Radial_b_pl(lambda,mu,nu)/anorm;
  elif radial_op=Radial_bm and lambda_var=1 then
    ME_Radial_bm_pl(lambda,mu,nu)*anorm;
  elif radial_op=Radial_Db and lambda_var=1 then
    ME_Radial_Db_pl(lambda,mu,nu)*anorm;

  elif radial_op=Radial_b and lambda_var=-1 then
    ME_Radial_b_ml(lambda,mu,nu)/anorm;
  elif radial_op=Radial_bm and lambda_var=-1 then
    ME_Radial_bm_ml(lambda,mu,nu)*anorm;
  elif radial_op=Radial_Db and lambda_var=-1 then
    ME_Radial_Db_ml(lambda,mu,nu)*anorm;

  elif radial_op=Radial_S0 and lambda_var=0 then
    ME_Radial_S0(lambda,mu,nu);
  elif radial_op=Radial_Sp and lambda_var=0 then
    ME_Radial_Sp(lambda,mu,nu);
  elif radial_op=Radial_Sm and lambda_var=0 then
    ME_Radial_Sm(lambda,mu,nu);


  elif radial_op=Radial_id and type(lambda_var,even) then
    if lambda_var>=0 then
      ME_Radial_id_pl(lambda,mu,nu,lambda_var/2):
    else
      ME_Radial_id_ml(lambda,mu,nu,-lambda_var/2):
    fi:

  else
      # might be a good idea to check that we have a valid radial operator
      if radial_op=Radial_id then
        MM:=RepRadial_Prod([],args[2..4],0,
                              max(mu,nu)+iquo(abs(lambda_var)+3,2)):
      else
        MM:=RepRadial_Prod([radial_op],args[2..4],0,
                              max(mu,nu)+iquo(abs(lambda_var)+3,2)):
      fi:
      MM[mu+1,nu+1]:  # matrices start at mu=nu=0!
  fi:

end:



# The following uses one of the above procedures to construct an explicit
# representation matrix for the corresponding lambda-changing operator.
# This matrix gives the action on the radial states labelled by
# nu_min..nu_max (and is thus of dimension nu_max-nu_min+1).
# It is the radial analogue of RepSO5_Y().
# The remember option is used, but these stored values are cleared
# each time the (much later) RepXspace_OpLC() routine is invoked.
# (Some of these matrices might be constructed more efficiently
#  somehow using a scan=diagonal option.) 

# Note that the datatype of the resulting matrix is not fixed:
# Maple chooses it - seemingly depending on the type of lambda
# (e.g. lambda=5/2 gives algebraic, lambda=2.5 gives floats;
#  these may be tested for using the Maple commands
#          type(MM,'Matrix'(datatype=anything));
#          type(MM,'Matrix'(datatype=float[8]));
# ).
# It'd thus be a good idea to apply evalf to all the matrix elements
# obtained before diagonalisation etc.

# FOR16: Typically, these routines will get called many times
#  during the construction of an operator (Hamiltonian) for
#  various values of lambda, but the same range of nu_min & nu_max.

RepRadial:=proc(ME::procedure,lambda::algebraic,
                                  nu_min::nonnegint,nu_max::nonnegint)
    option remember;

  Matrix(nu_max-nu_min+1,(i,j)->ME(lambda,nu_min-1+i,nu_min-1+j)):
end:

# The following returns the square root of the matrix obtained above.
# Again the remember option is used: and cleared by RepXspace_OpLC().

RepRadial_sq:=proc(ME::procedure,lambda::algebraic,
                                  nu_min::nonnegint,nu_max::nonnegint)
    option remember;

  MatrixFunction(evalf(RepRadial(ME,lambda,nu_min,nu_max)),sqrt(v),v):
end:


# Same as RepRadial above, but with an extra parameter to be passed.

RepRadial_param:=proc(ME::procedure,lambda::algebraic,
                           nu_min::nonnegint,nu_max::nonnegint,param::integer)
    option remember;

  Matrix(nu_max-nu_min+1,(i,j)->ME(lambda,nu_min-1+i,nu_min-1+j,param)):
end:


# The following takes a list of the lambda-changing radial operators,
# each one for which the change in beta is stipulated, and whose matrix
# elements are given by the above procedures.
# These operators are:
#    Radial_b2, Radial_bm2, Radial_D2b, Radial_bDb,  (lambda unchanged)
#    Radial_b_ml, Radial_bm_ml, Radial_Db_ml,  (lambda decremented)
#    Radial_b_sq, Radial_bm_sq, Radial_Db_sq,  (lambda unchanged)
#    Radial_b_pl, Radial_bm_pl, Radial_Db_pl,  (lambda incremented)
#    Radial_S0, Radial_Sp, Radial_Sm,  (lambda unchanged)
#    Radial_id   (lambda change according to parameter).
# The final operator here is the 1 operator, for which any even change
# in lambda is permitted. This change is specified by a parameter to
# the right of Radial_id in the list of operators. And to the right
# of that parameter, we expect to see the Radial_flag "operator" indicator.
# (We have organised things in this way to make parsing the operator
# list from the right easier).

# Thus lambda is carried along in the product (from right to left),
# so that the final value of lambda is determined by the list of operators.

RepRadialshift_Prod:=proc(rs_op::list(integer),
                             anorm::algebraic, lambda::algebraic,
                             nu_min::nonnegint, nu_max::nonnegint)
    local i,n,Mat,Mat_product,lambda_var,param;

  n:=nops(rs_op);

  if n=0 then    # null product: require identity matrix
    Mat_product:=Matrix([seq(1,i=nu_min..nu_max)],scan=diagonal):

  else
    lambda_var:=lambda:

    # form required product, multiplying from the right

    for i from n to 1 by -1 do

      # do inplace multiplications ...

      if rs_op[i]=Radial_b2 then
        Mat:=RepRadial(ME_Radial_b2,lambda_var,nu_min,nu_max);
        Mat:=MatrixScalarMultiply(Mat,1/anorm^2);
      elif rs_op[i]=Radial_bm2 then
        Mat:=RepRadial(ME_Radial_bm2,lambda_var,nu_min,nu_max);
        Mat:=MatrixScalarMultiply(Mat,anorm^2);
      elif rs_op[i]=Radial_D2b then
        Mat:=RepRadial(ME_Radial_D2b,lambda_var,nu_min,nu_max);
        Mat:=MatrixScalarMultiply(Mat,anorm^2);
      elif rs_op[i]=Radial_bDb then
        Mat:=RepRadial(ME_Radial_bDb,lambda_var,nu_min,nu_max);        

      elif rs_op[i]=Radial_b_pl then
        Mat:=RepRadial(ME_Radial_b_pl,lambda_var,nu_min,nu_max);
        lambda_var:=lambda_var+1:
        Mat:=MatrixScalarMultiply(Mat,1/anorm);
      elif rs_op[i]=Radial_bm_pl then
        Mat:=RepRadial(ME_Radial_bm_pl,lambda_var,nu_min,nu_max);
        lambda_var:=lambda_var+1:
        Mat:=MatrixScalarMultiply(Mat,anorm);
      elif rs_op[i]=Radial_Db_pl then
        Mat:=RepRadial(ME_Radial_Db_pl,lambda_var,nu_min,nu_max);
        lambda_var:=lambda_var+1:
        Mat:=MatrixScalarMultiply(Mat,anorm);

      elif rs_op[i]=Radial_b_ml then
        Mat:=RepRadial(ME_Radial_b_ml,lambda_var,nu_min,nu_max);
        lambda_var:=lambda_var-1:
        Mat:=MatrixScalarMultiply(Mat,1/anorm);
      elif rs_op[i]=Radial_bm_ml then
        Mat:=RepRadial(ME_Radial_bm_ml,lambda_var,nu_min,nu_max);
        lambda_var:=lambda_var-1:
        Mat:=MatrixScalarMultiply(Mat,anorm);
      elif rs_op[i]=Radial_Db_ml then
        Mat:=RepRadial(ME_Radial_Db_ml,lambda_var,nu_min,nu_max);
        lambda_var:=lambda_var-1:
        Mat:=MatrixScalarMultiply(Mat,anorm);

      elif rs_op[i]=Radial_b_sq then
          # obtain a matrix representing beta by taking the positive
          # definite square root of that representing beta^2.
        Mat:=RepRadial_sq(ME_Radial_b2,lambda_var,nu_min,nu_max):
        Mat:=MatrixScalarMultiply(Mat,1/anorm);
      elif rs_op[i]=Radial_bm_sq then
          # obtain a matrix representing 1/beta by taking the positive
          # definite square root of that representing 1/beta^2.
        Mat:=RepRadial_sq(ME_Radial_bm2,lambda_var,nu_min,nu_max):
        Mat:=MatrixScalarMultiply(Mat,anorm);
      elif rs_op[i]=Radial_Db_sq then
          # obtain a matrix representing d/d(beta) by taking the
          # positive definite square root of that representing 1/beta^2
          # multiplied by that for beta*d/d(beta).
        Mat:=MatrixMatrixMultiply(
                RepRadial_sq(ME_Radial_bm2,lambda_var,nu_min,nu_max),
                RepRadial(ME_Radial_bDb,lambda_var,nu_min,nu_max)):
        Mat:=MatrixScalarMultiply(Mat,anorm);

      elif rs_op[i]=Radial_S0 then
        Mat:=RepRadial(ME_Radial_S0,lambda_var,nu_min,nu_max);
      elif rs_op[i]=Radial_Sp then
        Mat:=RepRadial(ME_Radial_Sp,lambda_var,nu_min,nu_max);
      elif rs_op[i]=Radial_Sm then
        Mat:=RepRadial(ME_Radial_Sm,lambda_var,nu_min,nu_max);

      elif rs_op[i]=Radial_flag then
        if i<3 then
          error "operator flag %1 found at bad location %2", rs_op[i],i; 
        fi:

        param:=rs_op[i-1]:   # This should be even integer

        if rs_op[i-2]=Radial_id then
          if param>0 then
            Mat:=RepRadial_param(ME_Radial_id_pl,lambda_var,
                                     nu_min,nu_max,param/2);
          elif param<0 then
            Mat:=RepRadial_param(ME_Radial_id_ml,lambda_var,
                                     nu_min,nu_max,-param/2);
          fi:

          lambda_var:=lambda_var+param:

        else           # There's only one operator defined with parameter!
          error "operator %1 undefined with (flagged) parameter", rs_op[i-2]; 
        fi:

      else
        if rs_op[i]=Radial_id then
          error "operator %1 not used with (flagged) parameter", rs_op[i]; 
        else
          error "operator %1 undefined", rs_op[i]; 
        fi
      fi:

      if i=n then
        Mat_product:=Mat:
          # These matrices now have the same storage: but Mat_product is not
          # then changed when Mat is reassigned to another Matrix in the next
          # instance of loop.
      else
        Mat_product:=simplify(MatrixMatrixMultiply(Mat,Mat_product)):
      fi:

      if rs_op[i]=Radial_flag then   # need to skip two entries
        i:=i-2
      fi:
    od:
  
  fi:

  Mat_product:
end:

# Instead of specifying a list of lambda_changing operators, it will
# usually be the case that we wish to represent an arbitrary product
# of the "basic" radial operators
#   beta^2; 1/beta^2; d^2/d(beta)^2; beta*d/d(beta);
#   beta; 1/beta; d/d(beta); S0; S+; S-;
# between two spaces which have certain values of lambda
# (possibly different - but differing by an integer).

# The following procedure takes a list [a1,a2,a3,a4,.....], with each
# element one of the values from
#    Radial_b2, Radial_bm2, Radial_D2b, Radial_bDb,
#    Radial_b,  Radial_bm, Radial_Db,
#    Radial_S0, Radial_Sp, Radial_Sm,
# and lambda_var which is the difference in lambda between the two spaces.

# The list of basic radial operators is converted into a list of the
# lambda-changing operators, in such a way as to get the desired
# overall lambda change.
# We use the degree 1 operators as much as possible to effect a
# change in lambda, and split up a degree 2 operator
# (which would otherwise give no lambda change) if required.
# If the required lambda change is still not attained, then we
# insert a lambda-changing identity operator (which can only change
# lambda by an even number).
# If the parity of the resulting change is wrong, we use a non-analytic
# operator (which doesn't change lambda). In such a case, a matrix square
# root must be taken, and the resulting representation matrix but an
# apporixmation, which for a given number of
# states might only be reasonably accurate if far more states are used.

# (note that the su(1,1) operators are restricted here to keep lambda
#  unchanged: however, it ought to be possible to get a 2 change for them).

LambdaOp_List:=proc(rs_op::list(integer),lambda_var::integer)

  global Radial_List:
  local first_pass,i,opcount,opcount1,opcount2,split2,tosplit,splits,
        parity,lambda_run:

  # first get number of each of the radial operators

  opcount:=map2(numboccur,rs_op,Radial_List):
  if nops(rs_op)<>add(i,i in opcount) then
    error "undefined radial operator in %1", rs_op:
  fi:

  opcount2:=add(opcount[i],i=4..7):
  opcount1:=add(opcount[i],i=8..10):

#  if abs(lambda_var)>2*opcount2+opcount1 then
#    error "sweating to attain this lambda variation!":
#  fi:

  split2:=abs(lambda_var)-opcount1:
  if split2<=0 then
    first_pass:=rs_op:  # no splitting of deg 2 operators required.

  else
    tosplit:=iquo(split2+1,2):  # need to split up this num of deg2 operators

    splits:=[0,0,0,0]:  # num to split of each type

    splits[1]:=min(tosplit,opcount[4]):    # first split beta^2
    tosplit:=tosplit-splits[1]:
    splits[4]:=min(tosplit,opcount[7]):    # then split beta*d/d(beta)
    tosplit:=tosplit-splits[4]:
    splits[3]:=min(tosplit,opcount[6]):    # then split d^2/d(beta^2)
    tosplit:=tosplit-splits[3]:
    splits[2]:=min(tosplit,opcount[5]):    # then split beta^(-2)
    tosplit:=tosplit-splits[2]:

    #  create first_pass, splitting degree 2 ops as necessary

    first_pass:=[]:
    for i from 1 to nops(rs_op) do
      if rs_op[i] in [Radial_Sm,Radial_S0,Radial_Sp,
                      Radial_b, Radial_bm, Radial_Db] then
        first_pass:=[op(first_pass),rs_op[i]]:

      elif rs_op[i]=Radial_b2 then
        if splits[1]>0 then
          first_pass:=[op(first_pass),Radial_b,Radial_b]:
          splits[1]:=splits[1]-1:
        else
          first_pass:=[op(first_pass),Radial_b2]:
        fi:

      elif rs_op[i]=Radial_bDb then
        if splits[4]>0 then
          first_pass:=[op(first_pass),Radial_b,Radial_Db]:
          splits[4]:=splits[4]-1:
        else
          first_pass:=[op(first_pass),Radial_bDb]:
        fi:

      elif rs_op[i]=Radial_D2b then
        if splits[3]>0 then
          first_pass:=[op(first_pass),Radial_Db,Radial_Db]:
          splits[3]:=splits[3]-1:
        else
          first_pass:=[op(first_pass),Radial_D2b]:
        fi:

      elif rs_op[i]=Radial_bm2 then
        if splits[2]>0 then
          first_pass:=[op(first_pass),Radial_bm,Radial_bm]:
          splits[2]:=splits[2]-1:
        else
          first_pass:=[op(first_pass),Radial_bm2]:
        fi:
      fi:
    od:

  fi:

  # now make a pass through the operators, assigning a lambda direction
  # to the degree 1 operators.

  parity:=irem(split2,2); # if 1, we must use a degree 1 op with lambda const.

  # But if there are no operators, we must introduce a pair that gives
  # the identity operator (to eventually yield +/-1 lambda variation).

  if parity=1 and nops(first_pass)=0 then
    first_pass:=[Radial_b, Radial_bm]:
  fi:

  lambda_run:=lambda_var: # required lambda after (to left of) current operator

  for i to nops(first_pass) do
    if first_pass[i] in [Radial_b, Radial_bm, Radial_Db] then
      if parity<>0 then    # use constant lambda operator
        first_pass[i]:=first_pass[i]+Radial_shift_sq:
        parity:=0:
      elif lambda_run>0 then
        first_pass[i]:=first_pass[i]+Radial_shift_pl:
        lambda_run:=lambda_run-1:
      else # lambda_run<=0 then
        first_pass[i]:=first_pass[i]+Radial_shift_ml:
        lambda_run:=lambda_run+1:
      fi:
    fi:
  od:

  # now if we haven't yet reached the correct lambda, put in an
  # identity operator with the required variation
  # (best is at the front for +ve variation, at the back for -ve).

  if lambda_run>0 then
    [Radial_id,lambda_run,Radial_flag,op(first_pass)]:
  elif lambda_run<0 then
    [op(first_pass),Radial_id,lambda_run,Radial_flag]:
  else
    first_pass:
  fi:

end:




LambdaOp_List_Old:=proc(rs_op::list(integer),lambda_var::integer)

  global Radial_List:
  local first_pass,i,opcount,opcount1,opcount2,split2,tosplit,splits,
        parity,lambda_run:

  # first get number of each of the radial operators

  opcount:=map2(numboccur,rs_op,Radial_List):
  if nops(rs_op)<>add(i,i in opcount) then
    error "undefined radial operator in %1", rs_op:
  fi:

  opcount2:=add(opcount[i],i=4..7):
  opcount1:=add(opcount[i],i=8..10):

#  if abs(lambda_var)>2*opcount2+opcount1 then
#    error "this operator cannot attain this lambda variation":
#  fi:

  split2:=abs(lambda_var)-opcount1:
  if split2<=0 then
    first_pass:=rs_op:  # no splitting of deg 2 operators required.

  else
    tosplit:=iquo(split2+1,2):  # need to split up this num of deg2 operators

    splits:=[0,0,0,0]:  # num to split of each type

    if tosplit<=opcount[4] then    # first split beta^2
      splits[1]:=tosplit:
    else
      splits[1]:=opcount[4]:   
      tosplit:=tosplit-opcount[4]:

      if tosplit<=opcount[7] then    # then split beta*d/d(beta)
        splits[4]:=tosplit:
      else
        splits[4]:=opcount[7]:   
        tosplit:=tosplit-opcount[7]:

        if tosplit<=opcount[6] then    # then split d^2/d(beta^2)
          splits[3]:=tosplit:
        else
          splits[3]:=opcount[6]:   

          splits[2]:=tosplit-opcount[6]:  # rest in beta^(-2)
        fi:
      fi:
    fi:

    #  create first_pass, splitting degree 2 ops as necessary

    first_pass:=[]:
    for i from 1 to nops(rs_op) do
      if rs_op[i] in [Radial_Sm,Radial_S0,Radial_Sp,
                      Radial_b, Radial_bm, Radial_Db] then
        first_pass:=[op(first_pass),rs_op[i]]:

      elif rs_op[i]=Radial_b2 then
        if splits[1]>0 then
          first_pass:=[op(first_pass),Radial_b,Radial_b]:
          splits[1]:=splits[1]-1:
        else
          first_pass:=[op(first_pass),Radial_b2]:
        fi:

      elif rs_op[i]=Radial_bDb then
        if splits[4]>0 then
          first_pass:=[op(first_pass),Radial_b,Radial_Db]:
          splits[4]:=splits[4]-1:
        else
          first_pass:=[op(first_pass),Radial_bDb]:
        fi:

      elif rs_op[i]=Radial_D2b then
        if splits[3]>0 then
          first_pass:=[op(first_pass),Radial_Db,Radial_Db]:
          splits[3]:=splits[3]-1:
        else
          first_pass:=[op(first_pass),Radial_D2b]:
        fi:

      elif rs_op[i]=Radial_bm2 then
        if splits[2]>0 then
          first_pass:=[op(first_pass),Radial_bm,Radial_bm]:
          splits[2]:=splits[2]-1:
        else
          first_pass:=[op(first_pass),Radial_bm2]:
        fi:
      fi:
    od:

  fi:

  # now make a pass through the operators, assigning a lambda direction
  # to the degree 1 operators.

  parity:=irem(split2,2); # if 1, we must use a degree 1 op with lambda const.
  lambda_run:=lambda_var: # required lambda after (to left of) current operator

  for i to nops(first_pass) do
    if first_pass[i] in [Radial_b, Radial_bm, Radial_Db] then
      if parity<>0 then    # use constant lambda operator
        first_pass[i]:=first_pass[i]+Radial_shift_sq:
        parity:=0:
      elif lambda_run>0 then
        first_pass[i]:=first_pass[i]+Radial_shift_pl:
        lambda_run:=lambda_run-1:
      else # lambda_run<=0 then
        first_pass[i]:=first_pass[i]+Radial_shift_ml:
        lambda_run:=lambda_run+1:
      fi:
    fi:
  od:

  first_pass:

end:


# The following function has a new interface!!!!!!!!!

# The next function calculates the matrix representing an arbitrary
# product of the "basic" radial operators
#   beta^2; 1/beta^2; d^2/d(beta)^2; beta*d/d(beta);
#   beta; 1/beta; d/d(beta); S0; S+; S-;
# The matrix elements are between basis elements for which the
# change in lambda is as stipulated in the function arguments:
#               from lambda -> lambda+lambda_var.

# The result might need evalf operating on it to get a legible matrix result.

# This product is encoded as an array [a1,a2,a3,a4,.....], with each
# element one of the values from
#    Radial_b2, Radial_bm2, Radial_D2b, Radial_bDb,
#    Radial_b,  Radial_bm, Radial_Db,
#    Radial_S0, Radial_Sp, Radial_Sm.

# Implementation details:
# The procedure LambdaOp_List (above) chooses, as appropriate, how each
# of these operators is to change lambda, and obtains the particular
# representation matrix.
# Sometimes it will be necessary to use one of the odd degree operators
# with a 0 change in lambda. Then a matrix square root must be taken,
# and the result but an apporixmation, which for a given number of
# states might only be reasonably accurate if far more states are used.


RepRadial_Prod:=proc(rbs_op::list(integer), anorm::algebraic,
                          lambda::algebraic, lambda_var::integer,
                          nu_min::nonnegint, nu_max::nonnegint)
    local rs_op:

    rs_op:=LambdaOp_List(rbs_op,lambda_var):
    RepRadialshift_Prod(rs_op,anorm,lambda,nu_min,nu_max):
end;

# As the above function, but using remember table.

RepRadial_Prod_rem:=proc(rbs_op::list(integer), anorm::algebraic,
                          lambda::algebraic, lambda_var::integer,
                          nu_min::nonnegint, nu_max::nonnegint)
    option remember;
  RepRadial_Prod(args):
end;


# The next procedure is more general in that it provides linear
# combinations of such products.
# The arguments anorm, lambda, lambda_var, nu_min, nu_max are same as above,
# the first, rlc_op, is of the form
#         [ [coeff1,rs_op1], [coeff2,rs_op2], ...],
# where rs_op1, rs_op2 are lists of basic radial operators, as in the first
# argument above.
# This procedure is only used by the procedures that represent the
# intrinsic Xspace operators pi, pi x pi, pi x q x pi, i.e., by the
# procedures RepXspace_Pi, RepXspace_PiPi and RepXspace_PiqPi.


RepRadial_LC:=proc(rlc_op::list(list), anorm::algebraic,
                          lambda::algebraic, lambda_var::integer,
                          nu_min::nonnegint, nu_max::nonnegint)
    local i,n,Mat;

  n:=nops(rlc_op);

  if n=0 then
    Mat:=Matrix(nu_max-nu_min+1);  #Null matrix
  else
    Mat:=MatrixScalarMultiply(
            evalf(RepRadial_Prod_rem(rlc_op[1][2],anorm,
                                   lambda,lambda_var,nu_min,nu_max)),
            rlc_op[1][1]);
    for i from 2 to n do
      Mat:=MatrixAdd(Mat,
            evalf(RepRadial_Prod_rem(rlc_op[i][2],anorm,
                                   lambda,lambda_var,nu_min,nu_max)),
            1,rlc_op[i][1],inplace);
    od:
  fi:
    
  Mat;
end:


###########################################################################
####------------- Representing Operators on full Xspace ---------------####
###########################################################################

# The following procedure RepXspace_Prod() forms the (adjusted SO(3)-reduced)
# matrix representing a product of the above defined operators on the
# product space of the radial states nu_min..nu_max and the spherical
# states of seniorities v_min..v_max. (For operators that are linear
# combinations of products, see ProdLC_Mat below.)
# The correct 4*Pi factors are included.

# By adjusted, we mean that the matrix elements should be mulitplied by
# sqrt(2*L_f+1) to get the genuine SO(3)-reduced matrix elements of
# the operator in question. These are useful because in practical
# use, because this 1/sqrt(2*L_f+1) appears in the Wigner-Eckart theorem
# (as in eqn. (4.141) of [RowanWood]).
# In the case of Hamiltonians: then L=M=0 and (L_i M_i 0 0 | L_f M_f)=1
# and this return value returns the required amplitude directly. 

# Either a single angular momentum (L) is used, or a range (L..L_max).
# BEWARE with the L arguments: the L_max argument is optional - if omitted,
# the single value of L is used; and if L_max is given, all values between
# L and Lmax (inclusive) are used.
# The matrix elements are all floating point numbers.

## A lambda value is associated with each state.
## In this procedure, it depends only on the value of v mod 2.
## It is given by (lambda_base + v mod 2).
## (This is the case in the paper [RWC2009] - see section VII.B therein.
## In general, the value of lambda might vary more - but this is not
## dealt with here).
# The SO(5) tensor operators are not affected by this value of lambda.
# However, the radial operators are, and thus different versions of
# some radial operators are used depending on the value of v.
# In addition, odd powers of the beta operator will swap between the
# lambda values (this is reflected in the coding of the above procedure).

# x_ops is a list of integer values, each of which represents an operator,
# as described above.
# The list then denotes the product of these operators.
# The allowed operators are

# SpHarm_112     #  Y_{112}
# SpHarm_212     #  Y_{212}
# SpHarm_214     #  Y_{214}
# SpHarm_310     #  Y_{310}
# SpHarm_313     #  Y_{313}
# SpHarm_314     #  Y_{314}
# SpHarm_316     #  Y_{316}
# (SpHarm_vaL    #  for all v<=6)
# SpDiag_sqLdim  #  Diagonal matrix with entries (-1)^L*sqrt(2L+1).
# SpDiag_sqLdiv  #  Diagonal matrix with entries (-1)^L/sqrt(2L+1).
# Radial_Sm      #  SU(1,1) operator S-
# Radial_S0      #  SU(1,1) operator S0
# Radial_Sp      #  SU(1,1) operator S+
# Radial_b2      #  beta^2
# Radial_bm2     #  1/beta^2
# Radial_D2b     #  d^2/d(beta)^2
# Radial_bDb     #  beta*d/d(beta)
# Radial_b       #  beta
# Radial_bm      #  1/beta
# Radial_Db      #  d/d(beta)
# Xspace_Pi      #  For operator  pi;
# Xspace_PiPi2   #  For operator  [pi x pi]_{v=2,L=2};
# Xspace_PiPi4   #  For operator  [pi x pi]_{v=2,L=2};
# Xspace_PiqPi     #  For the "grubby" operator  [q x pi x pi]_{v=3,L=0};

# For an empty list, the identity matrix is returned.

# With regard to labelling the basis vectors of the tensor product space,
# the nu label varies quickest, and the sph labels [v,alpha,L] slowest.
# When L varies, it does so most slowest (so that the matrices formed
# for a range of L values are obtained by simply adjoining those obtained
# for the individual L values).
# This corresponds to the order of the state labels output by lbsXspace().


RepXspace_Prod:=proc(x_ops::list(integer),
                     anorm::algebraic,lambda_base::algebraic,
                     nu_min::nonnegint,nu_max::nonnegint,
                     v_min::nonnegint,v_max::nonnegint,
                     L::nonnegint,L_max::nonnegint,$)
    local sph_ops,nu_ops,run_Mat,xsp_Mat,this_op,up_running;
    global Radial_Max,Spherical_Min,
           Xspace_Pi,Xspace_PiPi2,Xspace_PiPi4,Xspace_PiqPi;

  # Run through the list of operators left to right, storing independently
  # those that act on the radial and spherical spaces.
  # If/When we see the grubby or pi or [pi x pi] operators,
  # we form the matrices and multiply them out.

  up_running:=0:            # flag to indicate that run_Mat contains sommat
  sph_ops:=[]: nu_ops:=[]:  # accumulates operators from the left
 
  for this_op in x_ops do

    if not type(this_op,integer) then
        error "Operator %1 undefined", this_op:
    elif this_op<=Radial_Max then
      nu_ops:=[op(nu_ops),this_op];     # store Radial Ops
    elif this_op>=Spherical_Min then
      sph_ops:=[op(sph_ops),this_op];   # store Sph Ops
    else

      # we now expect an operator on the full Xspace, so we need to multiply
      # out all those on the Radial and Sph spaces so far accumulated.

      if nu_ops<>[] or sph_ops<>[] then
        xsp_Mat:=RepXspace_Twin(nu_ops,sph_ops,args[2..-1]):
        nu_ops:=[]:   # used, so reset
        sph_ops:=[]:  # ditto
        if up_running>0 then
          MatrixMatrixMultiply(run_Mat,xsp_Mat,inplace):
        else          # nothing yet, so use xsp_Mat (not a copy)
          run_Mat:=xsp_Mat:
          up_running:=1:
        fi:
      fi:

      # generate the Xspace operator as required

      if this_op=Xspace_PiqPi then
          xsp_Mat:=RepXspace_PiqPi(args[2..-1]):

      elif this_op=Xspace_PiPi2 then   # For operator  [pi x pi]_{v=2,L=2};
          xsp_Mat:=RepXspace_PiPi(2,args[2..-1]):

      elif this_op=Xspace_PiPi4 then   # For operator  [pi x pi]_{v=2,L=4};
          xsp_Mat:=RepXspace_PiPi(4,args[2..-1]):

      elif this_op=Xspace_Pi then   # For operator  [pi]_{v=1,L=2};
          xsp_Mat:=RepXspace_Pi(args[2..-1]):

      # put other Xpsace operators here!

      else
        error "Operator %1 undefined", this_op:
      fi:

      # Now multiply in this Xspace operator

      if up_running>0 then
        MatrixMatrixMultiply(run_Mat,xsp_Mat,inplace);
      else
        run_Mat:=Matrix(xsp_Mat):   # need a copy because of remember tables
        up_running:=1:
      fi:
    fi:
  od:

  # And we must multiply out any remaining operators.

  if nu_ops<>[] or sph_ops<>[] then
    xsp_Mat:=RepXspace_Twin(nu_ops,sph_ops,args[2..-1]):

    if up_running>0 then
      MatrixMatrixMultiply(run_Mat,xsp_Mat,inplace):
    else          # nothing yet, so use xsp_Mat (not a copy)
      run_Mat:=xsp_Mat:
      up_running:=1:
    fi:
  fi:

  if up_running=0 then   # empty operator - need identity matrix
    run_Mat:=Matrix([seq(1,i=1..dimXspace(args[4..-1]))],scan=diagonal):
  fi:

  run_Mat:

end:


# The following function takes two lists of operators,
# given as two separate lists, the first of which only acts on
# the radial space (but with a dependence on seniority), and the
# second only acts on the spherical space.
# The matrix that is returned is the combined action on the
# cross-product space, taking account of the differing lambda values
# for the different seniorities.

# Note that the correct 4*Pi factors are applied here, so that the
# matrix elements returned are genuine adjusted SO(3)-reduced.


RepXspace_Twin:=proc(rad_ops::list(integer), sph_ops::list(integer),
                     anorm::algebraic, lambda_base::algebraic,
                     nu_min::nonnegint, nu_max::nonnegint,
                     v_min::nonnegint, v_max::nonnegint,
                     L::nonnegint, L_max::nonnegint,$)
    local j2,i2,j1,i1,jdisp,idisp,lambda_disp_init,lambda_disp_fin,
          rad_dim,rad_Mat,
          sph_dim,sph_labels,sph_Mat,sph_ME,
          direct_Mat;
    global glb_lam_fun;


  # Obtain dim, labels for S5 space, and the rep matrix of the sph operator

  sph_dim:=dimSO5r3_rngVvarL(args[7..-1]);
  sph_labels:=lbsSO5r3_rngVvarL(args[7..-1]);
  sph_Mat:=RepSO5r3_Prod(sph_ops,args[7..-1]);
  # Fix the (4pi) factors in the latter
  MatrixScalarMultiply(sph_Mat,
        evalf(Convert_red^NumSO5r3_Prod(sph_ops)),inplace);

  # dimension of radial space:

  rad_dim:=dimRadial(nu_min,nu_max);

  # now form the direct product representations on the space of the
  # following dimension.

  direct_Mat:=Matrix(sph_dim*rad_dim,datatype=float[8]); # setting element type

  # Place the entries one-by-one into the direct product matrix.
  # (j1;j2) is initial state, (i1;i2) is final state
  # 1st label is radial, 2nd is SO5, and varies slowest.

  # For the radial (nu) part, the rep matrix depends on the initial
  # and final values of v: these determine the lambdas mapped between.

  for j2 to sph_dim do
    jdisp:=(j2-1)*rad_dim:
    lambda_disp_init:=glb_lam_fun(sph_labels[j2][1]):

  for i2 to sph_dim do
    idisp:=(i2-1)*rad_dim:
    lambda_disp_fin:=glb_lam_fun(sph_labels[i2][1]):
    sph_ME:=sph_Mat[i2,j2]:

    if sph_ME=0 then
#      print("Omitting ",sph_labels[j2],":",sph_labels[i2]):
      next
    fi:  # skip zero cases of spherical MEs.

    rad_Mat:=RepRadial_Prod_rem(rad_ops,anorm,
                                  lambda_base+lambda_disp_init,
                                  lambda_disp_fin-lambda_disp_init,
                                  nu_min,nu_max):

    for i1 to rad_dim do
    for j1 to rad_dim do
      direct_Mat[idisp+i1,jdisp+j1]:=evalf(rad_Mat[i1,j1]*sph_ME);
    od: od:
  od: od:

  direct_Mat;
end:


# The following procedure obtains representation matrices in the
# product space as above, except that here, linear combinations of
# operators are possible. These are encoded in x_oplc.
# This parameter x_oplc is a list of the form [ [coeff1,x_ops1],
# [coeff2,x_ops2], ....].
# Here coeff1,coeff2 etc are coefficients for the operator products
# encoded in x_ops1, x_ops2 etc.
# Each x_ops is itself a list, comprising elements from the list
# given above RepXspace_Prod() above.

# If the list is empty, the null matrix is returned.

# The coefficients coeff1, coeff2, etc., can be numerical values OR they
# can be functions of SENIORITY, ANGMOM, NUMBER, ALFA (anything else
# will cause problems!) - these will be substituted for according to the
# [nu,v,alpha,L] values of the state being operated on by setting
# SENIORITY=v, ANGMOM=L, NUMBER=nu, ALFA=alpha
# (we probably will never use alfa!)

RepXspace_OpLC:=proc(x_oplc::list, anorm::algebraic, lambda_base::algebraic,
                     nu_min::nonnegint, nu_max::nonnegint,
                     v_min::nonnegint, v_max::nonnegint,
                     L::nonnegint, L_max::nonnegint,$)
      local Rmat,Pmat,i,n,Xlabels;

  # first clear the remember tables for the SO5 and radial matrices
  # so that we can start afresh with the given Xspace.

  forget(RepSO5_Y):
  forget(RepRadial_Prod_rem):
  forget(RepRadial):
  forget(RepRadial_sq):
  forget(RepRadial_param):
  forget(RepXspace_PiqPi):
  forget(RepXspace_PiPi):
  forget(RepXspace_Pi):
  n:=nops(x_oplc);

  Xlabels:=lbsXspace(args[4..-1]):    # list of all states in X-space.
 
  if n=0 then    # null sum: require zero matrix
    Rmat:=Matrix( dimXspace(args[4..-1]), datatype=float[8] ); #Null matrix
  else

    # first obtain rep matrix on X-space of 1st operator product

    Rmat:=RepXspace_Prod(x_oplc[1][2],args[2..-1]);

    if type(x_oplc[1][1],constant) then

      # simply multiply by the coefficient (which is a numeric value)

      MatrixScalarMultiply(Rmat,evalf(x_oplc[1][1]),inplace);
    else

      # post-multiply by a diagonal matrix that is formed by evaluating
      # the coefficient (a function of number, seniority, alfa, angmom)
      # at each state in the X_space

      MatrixMatrixMultiply(Rmat,
            Matrix(map(x->evalf(eval(x_oplc[1][1],
                          [NUMBER=x[1],SENIORITY=x[2],ALFA=x[3],ANGMOM=x[4]])),
                       Xlabels),
                       shape=diagonal,scan=diagonal), inplace);
    fi:

    # now do similar for every other operator product - and sum results

    for i from 2 to n do
      if type(x_oplc[i][1],constant) then
        MatrixAdd(Rmat, RepXspace_Prod(x_oplc[i][2],args[2..-1]),
                                     1, evalf(x_oplc[i][1]),inplace);
      else
        Pmat:=RepXspace_Prod(x_oplc[i][2],args[2..-1]):
        MatrixMatrixMultiply(Pmat,
              Matrix(map(x->evalf(eval(x_oplc[i][1],
                         [NUMBER=x[1],SENIORITY=x[2],ALFA=X[3],ANGMOM=x[4]])),
                         Xlabels),
                         shape=diagonal,scan=diagonal),inplace);
        MatrixAdd(Rmat,Pmat,inplace);
      fi:
    od:
  fi:

  Rmat;
end:



# The following provides a representation of the pi/(-i\hbar) operator.
# on the product space.

# If an exotic coefficient is required
# (e.g. a function of NUMBER, SENIORITY, ANGMOM or ALFA),
# the procedure RepXspace_OpLC can be used to call this one.

# The returned matrix elements are ADJUSTED
# (i.e. lack the sqrt(2L_f+1)) SO(3)-reduced matrix elements.

RepXspace_Pi:=proc( anorm::algebraic, lambda_base::algebraic,
                    nu_min::nonnegint, nu_max::nonnegint,
                    v_min::nonnegint, v_max::nonnegint,
                    L::nonnegint, L_max::nonnegint,$)
    option remember;
    local j2,i2,j1,i1,jdisp,idisp,lambda_disp_init,lambda_disp_fin,
          v_init,al_init,L_init,v_fin,al_fin,L_fin,v_chg,
          sph_dim,sph_labels,
          rad_dim,rad_Mat,direct_Mat,CG2;
    global glb_lam_fun,Radial_Db,Radial_bm;

  # dimension of radial space:

  rad_dim:=dimRadial(nu_min,nu_max);

  # Obtain dim and labels for S5 space.

  sph_dim:=dimSO5r3_rngVvarL(args[5..-1]);
  sph_labels:=lbsSO5r3_rngVvarL(args[5..-1]);

  # Will form the representation on the sph_dim*rat_dim dimensional
  # direct product space in the following Matrix, which is returned.
  # Using datatype=float[8] is necessary for later!

  direct_Mat:=Matrix(sph_dim*rad_dim,datatype=float[8]);

  # Place the entries one-by-one into the direct product matrix.
  # (j1;j2) is initial state, (i1;i2) is final state.
  # 1st label is radial, 2nd is SO5, and varies slowest.

  # For the radial (nu) part, the rep matrix depends on the initial
  # and final values of v: these determine the lambdas mapped between.

  for j2 to sph_dim do
    v_init:=sph_labels[j2][1]:   # seniority of initial state
    al_init:=sph_labels[j2][2]:  # alpha of same
    L_init:=sph_labels[j2][3]:   # L of same
    jdisp:=(j2-1)*rad_dim:
    lambda_disp_init:=glb_lam_fun(v_init):

  for i2 to sph_dim do
    L_fin:=sph_labels[i2][3]:    # L of final state
    if L_fin-L_init>2 or L_fin-L_init<-2 then next fi:  # zero becos q has L=2.

    v_fin:=sph_labels[i2][1]:    # seniority of final state
    al_fin:=sph_labels[i2][2]:   # alpha of same
    v_chg:=v_fin-v_init:         # change in seniority
    idisp:=(i2-1)*rad_dim:
    lambda_disp_fin:=glb_lam_fun(v_fin):


    # Now obtain the (SO(5) reduced) representation matrix between these
    # subspaces having constant spherical labels by treating separately
    # the cases v_chg = 1 (using (47a)) and -1 (using (47b)).

    if v_chg = 1 then

       rad_Mat:=RepRadial_LC([ [1,[Radial_Db]],[-v_init-2,[Radial_bm]] ],
                   anorm, lambda_base+lambda_disp_init,
                   lambda_disp_fin-lambda_disp_init, nu_min, nu_max):

      # multiply this rad_Mat radial action by the SO(5) reduced ME
      # from (18).

      MatrixScalarMultiply(rad_Mat,evalf(Qred_p1(v_init)),inplace);

    elif v_chg = -1 then

       rad_Mat:=RepRadial_LC([ [1,[Radial_Db]],[v_init+1,[Radial_bm]] ],
                   anorm, lambda_base+lambda_disp_init,
                   lambda_disp_fin-lambda_disp_init, nu_min, nu_max):

      # multiply this nu_Mat radial action by the SO(5) reduced ME

      MatrixScalarMultiply(rad_Mat,evalf(Qred_m1(v_init)),inplace);

    else
      next           # skip this (j1,i1) case (because it is zero).
    fi:

    # The "adjusted" SO(3)-reduced matrix element is obtained from the
    # SO(5)-reduced ME calculated above by multiplying with the following
    # (we could multiply by sqrt(dimSO3(L_fin)) here to get genuine
    #  SO(3)-reduced matrix elements):

    CG2:=CG_SO5r3(v_init,al_init,L_init,1,1,2,v_fin,al_fin,L_fin):

    # fill in the Xspace elements for these constant spherical
    # parameters (i2,j2).

    for i1 to rad_dim do
    for j1 to rad_dim do
      direct_Mat[idisp+i1,jdisp+j1]:=evalf(CG2*rad_Mat[i1,j1]);
    od: od:
  od: od:

  direct_Mat:
end:


# The following provides a representation of the operators
#            [pi x pi]_(v=2,L=2)   and  [pi x pi]_(v=2,L=4)
#            -------------------        -------------------
#                  hbar^2                    hbar^2
# on the product space (for these, PiPi_L is then 2 or 4).
# (I now believe I have the sign correct here!)

# If a more exotic coefficient is required
# (e.g. a function of NUMBER, SENIORITY, ANGMOM or ALFA),
# the procedure RepXspace_OpLC can be used to call this one.

# The returned matrix elements are ADJUSTED
# (i.e. lack the sqrt(2L_f+1)) SO(3)-reduced matrix elements.


RepXspace_PiPi:=proc(PiPi_L::nonnegint,
                     anorm::algebraic, lambda_base::algebraic,
                     nu_min::nonnegint, nu_max::nonnegint,
                     v_min::nonnegint, v_max::nonnegint,
                     L::nonnegint, L_max::nonnegint,$)
    option remember;
    local j2,i2,j1,i1,jdisp,idisp,lambda_disp_init,lambda_disp_fin,
          v_init,al_init,L_init,v_fin,al_fin,L_fin,v_chg,
          sph_dim,sph_labels,
          rad_dim,rad_Mat,direct_Mat,CG2;
    global glb_lam_fun,Radial_D2b,Radial_bm2,Radial_bDb;

  # dimension of radial space:

  rad_dim:=dimRadial(nu_min,nu_max);

  # Obtain dim and labels for S5 space.

  sph_dim:=dimSO5r3_rngVvarL(args[6..-1]);
  sph_labels:=lbsSO5r3_rngVvarL(args[6..-1]);

  # Will form the representation on the sph_dim*rat_dim dimensional
  # direct product space in the following Matrix, which is returned.
  # Using datatype=float[8] is necessary for later!

  direct_Mat:=Matrix(sph_dim*rad_dim,datatype=float[8]);

  # Place the entries one-by-one into the direct product matrix.
  # (j1;j2) is initial state, (i1;i2) is final state.
  # 1st label is radial, 2nd is SO5, and varies slowest.

  # For the radial (nu) part, the rep matrix depends on the initial
  # and final values of v: these determine the lambdas mapped between.

  for j2 to sph_dim do
    v_init:=sph_labels[j2][1]:   # seniority of initial state
    al_init:=sph_labels[j2][2]:  # alpha of same
    L_init:=sph_labels[j2][3]:   # L of same
    jdisp:=(j2-1)*rad_dim:
    lambda_disp_init:=glb_lam_fun(v_init):

  for i2 to sph_dim do
    L_fin:=sph_labels[i2][3]:    # L of final state
    if L_fin-L_init>PiPi_L or L_fin-L_init<-PiPi_L then next fi: # zero

    v_fin:=sph_labels[i2][1]:    # seniority of final state
    al_fin:=sph_labels[i2][2]:   # alpha of same
    v_chg:=v_fin-v_init:         # change in seniority
    idisp:=(i2-1)*rad_dim:
    lambda_disp_fin:=glb_lam_fun(v_fin):


    # Now obtain the (SO(5) reduced) representation matrix on this
    # subspace with constant spherical labels by treating separately
    # the cases v_chg = 2 (using (A10)) and -2 (using (A11))
    # and 0 (using (A12)).

    if v_chg = 2 then

       rad_Mat:=RepRadial_LC( [ [1,[Radial_D2b]],
                                   [(v_init+2)*(v_init+4),[Radial_bm2]],
                                   [-2*v_init-5,[Radial_bm2,Radial_bDb]] ],
                               anorm, lambda_base+lambda_disp_init,
                               lambda_disp_fin-lambda_disp_init,
                               nu_min, nu_max):

      # multiply this nu_Mat radial action by the SO(5) reduced ME
      # from (A2) (the minus sign comes from the i^2 (*hbar^2) )

      MatrixScalarMultiply(rad_Mat,-evalf(QxQred_p2(v_init)),inplace);

    elif v_chg = -2 then

       rad_Mat:=RepRadial_LC( [ [1,[Radial_D2b]],
                                   [(v_init-1)*(v_init+1),[Radial_bm2]],
                                   [2*v_init+1,[Radial_bm2,Radial_bDb]] ],
                               anorm, lambda_base+lambda_disp_init,
                               lambda_disp_fin-lambda_disp_init,
                               nu_min, nu_max):

      # multiply this nu_Mat radial action by the SO(5) reduced ME

      MatrixScalarMultiply(rad_Mat,-evalf(QxQred_m2(v_init)),inplace);

    elif v_chg = 0 then

       rad_Mat:=RepRadial_LC( [ [1,[Radial_D2b]],
                                   [-(v_init+1)*(v_init+2),[Radial_bm2]] ],
                               anorm, lambda_base+lambda_disp_init,
                               lambda_disp_fin-lambda_disp_init,
                               nu_min, nu_max):

      # multiply this nu_Mat radial action by the SO(5) reduced ME

      MatrixScalarMultiply(rad_Mat,-evalf(QxQred_0(v_init)),inplace);
    else
      next           # skip this (j1,i1) case (because it is zero).
    fi:


    # The "adjusted" SO(3)-reduced matrix element is obtained from the
    # SO(5)-reduced ME calculated above by multiplying with the following
    # (we could multiply by sqrt(dimSO3(L_fin)) here to get genuine
    #  SO(3)-reduced matrix elements)
    # (PiPi_L should be 2 or 4 here, else we'll get 0):

    CG2:=CG_SO5r3(v_init,al_init,L_init,2,1,PiPi_L,v_fin,al_fin,L_fin):

    # fill in the Xspace elements for these constant spherical
    # parameters (i2,j2).

    for i1 to rad_dim do
    for j1 to rad_dim do
      direct_Mat[idisp+i1,jdisp+j1]:=evalf(CG2*rad_Mat[i1,j1]);
    od: od:
  od: od:

  # This now correct for PiPi_L=4, but the sign needs to change for PiPi_L=2.

  if PiPi_L=2 then
      MatrixScalarMultiply(direct_Mat,-1,inplace);
  fi:  

  direct_Mat:
end:



# The following provides a representation of the "Grubby" operator
#                            [q x pi x pi]_(v=3)
#                            -------------------
#                                 hbar^2
# on the product space. (I'm not sure I have the sign correct here!)
# The argument parity will usually be 0. This indicates that the even
# seniority states have lambda=lambda_norm, with odd states having
# lambda=lambda_norm+1. For parity=1, it is the other way around
# (this isn't described in DJR's notes - and may be superfluous!)
# Slightly modified from DJR's expressions.

# If a more exotic coefficient is required (e.g. a function of
# NUMBER, SENIORITY, ANGMOM or ALFA), the procedure
# RepXspace_OpLC can be used to call this one.

# The returned matrix elements are ADJUSTED
# (i.e. lack the sqrt(2L_f+1)) SO(3)-reduced matrix elements.
# The result deviates from being hermitian due to boundary effects.


RepXspace_PiqPi:=proc( anorm::algebraic, lambda_base::algebraic,
                       nu_min::nonnegint, nu_max::nonnegint,
                       v_min::nonnegint, v_max::nonnegint,
                       L::nonnegint, L_max::nonnegint,$)
    option remember;
    local j2,i2,j1,i1,jdisp,idisp,lambda_disp_init,lambda_disp_fin,
          v_init,al_init,L_init,v_fin,al_fin,L_fin,v_chg,
          sph_dim,sph_labels,
          rad_dim,rad_Mat,rad_Mat2,direct_Mat,CG2;
    global glb_lam_fun,Radial_bm2,Radial_D2b,Radial_bDb,
                       Radial_b,Radial_Db,Radial_bm;

  # dimension of radial space:

  rad_dim:=dimRadial(nu_min,nu_max);

  # Obtain dim and labels for S5 space.

  sph_dim:=dimSO5r3_rngVvarL(args[5..-1]);
  sph_labels:=lbsSO5r3_rngVvarL(args[5..-1]);

  # Will form the representation on the sph_dim*rat_dim dimensional
  # direct product space in the following Matrix, which is returned.
  # Using datatype=float[8] is necessary for later!

  direct_Mat:=Matrix(sph_dim*rad_dim,datatype=float[8]);

  # Place the entries one-by-one into the direct product matrix.
  # (j1;j2) is initial state, (i1;i2) is final state.
  # 1st label is radial, 2nd is SO5, and varies slowest.

  # For the radial (nu) part, the rep matrix depends on the initial
  # and final values of v: these determine the lambdas mapped between.

  for j2 to sph_dim do
    v_init:=sph_labels[j2][1]:   # seniority of initial state
    al_init:=sph_labels[j2][2]:  # alpha of same
    L_init:=sph_labels[j2][3]:   # L of same
    jdisp:=(j2-1)*rad_dim:
    lambda_disp_init:=glb_lam_fun(v_init):

  for i2 to sph_dim do
    L_fin:=sph_labels[i2][3]:    # L of final state
    if L_fin<>L_init then next fi:   # need L's equal for this operator
    v_fin:=sph_labels[i2][1]:    # seniority of final state
    al_fin:=sph_labels[i2][2]:   # alpha of same
    v_chg:=v_fin-v_init:         # change in seniority
    idisp:=(i2-1)*rad_dim:
    lambda_disp_fin:=glb_lam_fun(v_fin):

    # Now obtain the (SO(5) reduced) representation matrix on this
    # subspace with constant spherical labels by treating separately
    # the cases v_chg = 3 (using (A14)) and -3 (using (A15))
    # and 1 (using (A17)) and -1 (using (A18)).
    # The SO(5)>SO(3) CG coefficient is left until the end.

    if v_chg = 3 then

       rad_Mat:=RepRadial_LC( [ [1,[Radial_b,Radial_D2b]],
                                   [(v_init+2)*(v_init+4),[Radial_bm]],
                                   [-2*v_init-5,[Radial_Db]] ],
                               anorm, lambda_base+lambda_disp_init,
                               lambda_disp_fin-lambda_disp_init,
                               nu_min, nu_max):

      # multiply this rad_Mat radial action by the SO(5) reduced ME
      # from (A4) (the minus sign in (A4) cancels that in (A14), which
      # comes from the i^2 (*hbar^2) ).

      MatrixScalarMultiply(rad_Mat,evalf(QxQxQred_p3(v_init)),inplace);

    elif v_chg = -3 then

       rad_Mat:=RepRadial_LC( [ [1,[Radial_b,Radial_D2b]],
                                   [(v_init-1)*(v_init+1),[Radial_bm]],
                                   [2*v_init+1,[Radial_Db]] ],
                               anorm, lambda_base+lambda_disp_init,
                               lambda_disp_fin-lambda_disp_init,
                               nu_min, nu_max):

      # multiply this rad_Mat radial action by the SO(5) reduced ME

      MatrixScalarMultiply(rad_Mat,evalf(QxQxQred_m3(v_init)),inplace);

    elif v_chg = 1 then

       rad_Mat:=RepRadial_LC( [ [1,[Radial_b,Radial_D2b]],
                                   [-(v_init+1)*(v_init+2),[Radial_bm]] ],
                               anorm, lambda_base+lambda_disp_init,
                               lambda_disp_fin-lambda_disp_init,
                               nu_min, nu_max):
       rad_Mat2:=RepRadial_LC( [ [1,[Radial_Db]],
                                 [-v_init-2,[Radial_bm]] ],
                               anorm, lambda_base+lambda_disp_init,
                               lambda_disp_fin-lambda_disp_init,
                               nu_min, nu_max):

      # combine rad_Mat and rad_Mat2 radials with appropriate SO(5) reduced MEs
      # (need to check signs!)

      MatrixAdd(rad_Mat,rad_Mat2,
                 +evalf(QxQxQred_p1(v_init)),
                 -evalf((2*v_init+5)*QixQxQred(v_init,v_fin,v_fin+1)),
                 inplace);

    elif v_chg = -1 then

       rad_Mat:=RepRadial_LC( [ [1,[Radial_b,Radial_D2b]],
                                   [-(v_init+1)*(v_init+2),[Radial_bm]] ],
                               anorm, lambda_base+lambda_disp_init,
                               lambda_disp_fin-lambda_disp_init,
                               nu_min, nu_max):
       rad_Mat2:=RepRadial_LC( [ [1,[Radial_Db]],
                                 [v_init+1,[Radial_bm]] ],
                               anorm, lambda_base+lambda_disp_init,
                               lambda_disp_fin-lambda_disp_init,
                               nu_min, nu_max):

      # combine rad_Mat and rad_Mat2 radials with appropriate SO(5) reduced MEs
      # (need to check signs!)

      MatrixAdd(rad_Mat,rad_Mat2,
                 +evalf(QxQxQred_m1(v_init)),
                 +evalf((2*v_init+1)*QixQxQred(v_init,v_fin,v_fin-1)),
                 inplace);

    else
      next           # skip this (j1,i1) case (because it is zero).
    fi:

    # The "adjusted" SO(3)-reduced matrix element is obtained from the
    # SO(5)-reduced ME calculated above by multiplying with the following
    # (we could multiply by sqrt(dimSO3(L_fin)) here to get genuine
    #  SO(3)-reduced matrix elements)

    CG2:=CG_SO5r3(v_init,al_init,L_init,3,1,0,v_fin,al_fin,L_fin):

    # fill up the Xspace elements for these constant spherical
    # parameters (i2,j2).

    for i1 to rad_dim do
    for j1 to rad_dim do
      direct_Mat[idisp+i1,jdisp+j1]:=evalf(CG2*rad_Mat[i1,j1]);
    od: od:
  od: od:

  direct_Mat:
end:



# The following procedure calculates genuine SO(5) reduced matrix elements
# for [Q^+ x Q x Q]_{v=3} and [Q^- x Q x Q]_{v=3}, which change v by +-1.
# This gives the (-)RHS of (A20), by evaluating a special case of the LHS.
# All four possible sign combinations work, although only two are
# needed for [pi x q x pi]_0. v_int is the "intermediate" v' value in (A20)
# Summing over all possible v_int should give -<v_f|||QxQxQ|||v_i>.
# In effect, this calculates (A23) & (A24) divided by (-)sqrt(2L+1) and
# the SO(5)>SO(3) CG coefficient (v,a_i,L;310||v_f,a_f,L).

QixQxQred:=proc(v_i::nonnegint,v_f::nonnegint,v_int::nonnegint)
  local mediates,L_i;

  L_i:=2*min(v_i,v_f):  # initial and final value of L

  # obtain the list of intermediate states with seniority v_int

  mediates:=lbsSO5r3_rngL(v_int,L_i-2,L_i+2):

  # Using (A19), sum over them to get a singly reduced matrix element.
  # Note that Q=4*Pi/sqrt(15) * Y112; [QxQ]_2=-4*Pi*sqrt(2/105) * Y212.

  ME_SO5red(v_f,1,v_int) * ME_SO5red(v_int,2,v_i) * sqrt(2/7) / 15
    * add( CG_SO5r3(v_int,m[2],m[3],1,1,2,v_f,1,L_i)
          * CG_SO5r3(v_i,1,L_i,2,1,2,v_int,m[2],m[3])
          * sqrt(dimSO3(m[3]))  * (-1)^m[3], m in mediates)

  # then convert this to a doubly reduced matrix element by dividing

    /CG_SO5r3(v_i,1,L_i,3,1,0,v_f,1,L_i)/sqrt(5*dimSO3(L_i)):

end;


# Previously the above four cases were implemented separately,
# by the following procedures. The correspondence is
#    QpxQxQred_p1(v)  <->  QixQxQred(v,v+1,v)
#    QmxQxQred_p1(v)  <->  QixQxQred(v,v+1,v+2)
#    QpxQxQred_m1(v)  <->  QixQxQred(v,v-1,v-2)
#    QmxQxQred_m1(v)  <->  QixQxQred(v,v-1,v)

# Note that, although all v>=0 are accepted as arguments for each,
# the first two don't make physical sense for v=0,
# and the last two don't make physical sense for v<=1 (both give errors).
# However, these exceptional values won't be required by Gruby_Op.

QpxQxQred_p1:=proc(v::nonnegint)   # via (A20) & (A22)
  local mediates,L1;

  L1:=2*v:  # initial value of L

  # obtain the list of intermediate states

  mediates:=lbsSO5r3_rngL(v,L1-2,L1):

  # sum over them to get a singly reduced matrix element

  Qred_p1(v) * QxQred_0(v) *
    add( CG_SO5r3(v,m[2],m[3],1,1,2,v+1,1,L1)
          * CG_SO5r3(v,1,L1,2,1,2,v,m[2],m[3])
          * sqrt(dimSO3(m[3]))  * (-1)^m[3], m in mediates)

  # then convert this to a doubly reduced matrix element by dividing

    /CG_SO5r3(v,1,L1,3,1,0,v+1,1,L1)/sqrt(5*dimSO3(L1)):

end;


QmxQxQred_p1:=proc(v::nonnegint)   # via (A20) & (A23)
  local mediates,L1;

  L1:=2*v:  # initial value of L

  # obtain the list of intermediate states

  mediates:=lbsSO5r3_rngL(v+2,L1-2,L1+2):

  # sum over them to get a singly reduced matrix element

  Qred_m1(v+2) * QxQred_p2(v) *
    add( CG_SO5r3(v+2,m[2],m[3],1,1,2,v+1,1,L1)
          * CG_SO5r3(v,1,L1,2,1,2,v+2,m[2],m[3])
          * sqrt(dimSO3(m[3])) * (-1)^m[3] , m in mediates)

  # then convert this to a doubly reduced matrix element by dividing

    /CG_SO5r3(v,1,L1,3,1,0,v+1,1,L1)/sqrt(5*dimSO3(L1)):

end;


QpxQxQred_m1:=proc(v::posint)   # via (77) & (80)
  local mediates,L1;

  L1:=2*v-2:  # initial value of L

  # this case has only one intermediate state [v-1,1,2v-4]

  # get singly reduced matrix element

  Qred_p1(v-2) * QxQred_m2(v)
          * CG_SO5r3(v-2,1,L1-2,1,1,2,v-1,1,L1)
          * CG_SO5r3(v,1,L1,2,1,2,v-2,1,L1-2)
          * sqrt(dimSO3(L1-2))

  # then convert this to a doubly reduced matrix element by dividing

    /CG_SO5r3(v,1,L1,3,1,0,v-1,1,L1)/sqrt(5*dimSO3(L1)):

end;


QmxQxQred_m1:=proc(v::posint)   # via (77) & (81)
  local mediates,L1;

  L1:=2*v-2:  # initial value of L

  # obtain the list of intermediate states

  mediates:=lbsSO5r3_rngL(v,L1-2,L1+2):

  # sum over them to get a singly reduced matrix element

  Qred_m1(v) * QxQred_0(v) *
    add( CG_SO5r3(v,m[2],m[3],1,1,2,v-1,1,L1)
          * CG_SO5r3(v,1,L1,2,1,2,v,m[2],m[3])
          * sqrt(dimSO3(m[3])) * (-1)^m[3] , m in mediates)

  # then convert this to a doubly reduced matrix element by dividing

    /CG_SO5r3(v,1,L1,3,1,0,v-1,1,L1)/sqrt(5*dimSO3(L1)):

end;




###########################################################################
####------------------- Specification of Operators --------------------####
###########################################################################


# The x_oplc parameter in the procedure RepXspace_OpLC() above
# depicts an operator on the R^5 space, that is composed from the
# operators on the radial (beta) and angular (S^5) spaces that
# have been described above. We may specify this parameter by hand,
# making use of the symbols defined for the various operators.
# However, it is simpler to use those defined next, or to use
# the procedure Hamiltonian_OpLC() further below.

# Here, we specify the Laplacian operator,
# and a couple of combinations of operators that
# may be used as a check on commutation relations:
# comm_su11_op defines a sum of three SU(1,1) operators
# which should produce a 0 matrix; comm_bdb_op similarly checks
# [d/d(beta) * beta, beta * d/d(beta)] respectively (note that an
# empty operator product [] denotes the identity operator).


laplacian_op:=[ [1,[Radial_D2b]],
                     [-(2+SENIORITY*(SENIORITY+3)),[Radial_bm2]] ]:


comm_su11_op:=[ [ 1,[Radial_Sm,Radial_Sp]],
                     [-1,[Radial_Sp,Radial_Sm]],
                     [-2,[Radial_S0]] ]:

comm_bdb_op:=[ [-1,[Radial_b,Radial_Db]],
                    [1,[Radial_Db,Radial_b]],
                    [-1,[]] ]:


# The hamiltonian_OpLC() (having up to 14 parameters) procedure produces
# the x_oplc operator for use in the procedure RepXspace_OpLC() above.
# This automatically creates lots of particular kinds of Hamitonians,
# so that the user does not have to bother with this internal code.

# In order, the up to 14 arguments specify coefficients of
#    Laplacian, 1, B^2, B^4,B^(-2),
#      B*cos(3g), B^3*cos(3g), B^5*cos(3g), B^(-1)*cos(3g),
#      cos(3g)^2, B^2*cos(3g)^2, B^4*cos(3g)^2, B^(-2)*cos(3g)^2,
#      [q x pi x pi]_(v=3).
# (Note that this ordering differs from a previous version!)

# In addition to simple numerical constants, these coefficients
# can be functions of SENIORITY, ANGMOM, NUMBER, ALFA
# (anything else will cause problems!) -
# these will be substituted for according to the [nu;v,alpha,L] values
# of the state being operated on by setting SENIORITY=v, ANGMOM=L,
# NUMBER=nu, ALFA=alpha.

Hamiltonian_OpLC:=proc(c11,c20,c21,c22,c23,
                           c30,c31,c32,c33,
                           c40,c41,c42,c43,
                           c50,$)
    local our_op:

  if nargs>0 and c11<>0 then   # build laplacian using eqn (69) of [RWC2009]
    our_op:=[ [c11,[Radial_D2b]],
              [-c11*(2+SENIORITY*(SENIORITY+3)),[Radial_bm2]] ]:
  else
    our_op:=[]:
  fi:

  if nargs>1 and c20<>0 then our_op:=[ op(our_op),
              [c20,[]] ]: fi:
  if nargs>2 and c21<>0 then our_op:=[ op(our_op),
              [c21,[Radial_b2]] ]: fi:
  if nargs>3 and c22<>0 then our_op:=[ op(our_op),
              [c22,[Radial_b2,Radial_b2]] ]: fi:
  if nargs>4 and c23<>0 then our_op:=[ op(our_op),
              [c23,[Radial_bm2]] ]: fi:
  if nargs>5 and c30<>0 then our_op:=[ op(our_op),
              [c30*Convert_310,[Radial_b,SpHarm_310]] ]: fi:
  if nargs>6 and c31<>0 then our_op:=[ op(our_op),
              [c31*Convert_310,[Radial_b2,Radial_b,SpHarm_310]] ]: fi:
  if nargs>7 and c32<>0 then our_op:=[ op(our_op),
              [c32*Convert_310,
                   [Radial_b2,Radial_b2,Radial_b,SpHarm_310]] ]: fi:
  if nargs>8 and c33<>0 then our_op:=[ op(our_op),
              [c33*Convert_310,[Radial_bm2,Radial_b,SpHarm_310]] ]: fi:
  if nargs>9 and c40<>0 then our_op:=[ op(our_op),
              [c40*Convert_310^2,[SpHarm_310,SpHarm_310]] ]: fi:
  if nargs>10 and c41<>0 then our_op:=[ op(our_op),
              [c41*Convert_310^2,[Radial_b2,SpHarm_310,SpHarm_310]] ]: fi:
  if nargs>11 and c42<>0 then our_op:=[ op(our_op),
              [c42*Convert_310^2,
                 [Radial_b2,Radial_b2,SpHarm_310,SpHarm_310]] ]: fi:
  if nargs>12 and c43<>0 then our_op:=[ op(our_op),
              [c43*Convert_310^2,[Radial_bm2,SpHarm_310,SpHarm_310]] ]: fi:
  if nargs>13 and c50<>0 then our_op:=[ op(our_op),
              [c50,[Xspace_PiqPi]] ]: fi:

  our_op:
end:


###########################################################################
####---------- Diagonalisation and Eigenbasis transformation ----------####
###########################################################################

# The following procedure DigXspace_OpLC() makes use of Eigenfiddle below
# to diagonalise the operator (Hamiltonian) ham_op, doing so independently
# on each L space, recording the eigenenergies and transformation matrix
# in each case.

# The return is a quartet of values
#           [eigen_vals, eigen_bases, Xparams, goodL];
# Here Xparams lists the parameters [lambda,anorm,nu_min,nu_max,v_min,v_max]
# without L. Here, goodL is a list of the values of angular momentum L
# between L_min..L_max which give non-zero dimensional spaces.
# The elements of the other two values pertain to these values of L.
# eigen_vals is a list each element of which contains a list of
# eigenvalues in a constant L space.
# eigen_bases is the list of transformation matrices to the eigenspaces.
# It's probably not a good idea to display this output!

# The values returned in eigen_vals are probably best displayed using
# the Show_Eigs() procedure below.
# The other output values may be directly used as input to the
# procedure AmpXspeig_OpLC() to calculate transition rates/amplitudes.

DigXspace_OpLC:=proc(ham_op::list,
                     anorm::algebraic, lambda_base::algebraic,
                     nu_min::nonnegint, nu_max::nonnegint,
                     v_min::nonnegint, v_max::nonnegint,
                     L_min::nonnegint, L_max::nonnegint,$)
      local LL,LLM,eigen_bases,eigen_result,goodL,
            L_count,Xparams,eigen_vals,h_matrix;

  if nargs=8 then  # no L_max
    LLM:=L_min;
  else
    LLM:=L_max
  fi:

  Xparams:=[anorm,lambda_base,nu_min,nu_max,v_min,v_max];  # 6 params (no L)

  # use goodL to store those L with non-zero dimension.
  # then only diagonalise these spaces.

  goodL:=[]:
  eigen_vals:=[];       # for eigenvalues for these L
  eigen_bases:=[];      # corresponding matrices of column vectors

  for LL from L_min to LLM do
    if dimSO5r3_rngV(v_min,v_max,LL)>0 then
      goodL:=[op(goodL),LL];

      h_matrix:=RepXspace_OpLC(ham_op,op(Xparams),LL);
      eigen_result:=Eigenfiddle(h_matrix);

      # store eigenvalue lists, one LL at a time.

      eigen_vals:=[op(eigen_vals),eigen_result[1]];

      # store matrix of eigenvectors (could do inverses here?)

      eigen_bases:=[op(eigen_bases),eigen_result[2]];
    fi:
  od:

  [eigen_vals, eigen_bases, Xparams, goodL];
end;



# The matrix H that is input to Eigenfiddle is diagonalised
# (using Maple's Eigenvectors procedure).
# This matrix is expected to be real.
# Before being diagonalised, it is averaged with its transpose
# to ensure that it is symmetric.
# The eigenvalues are returned in increasing order.

# The returned value is a pair:
#   [eigs_list,basis_Mat]
# the first element of this pair is a list of the lowest real eigenvalues;
# the second element of the pair is a matrix P which transforms the
# original matrix H to a diagonal matrix P^{-1}HP whose diagonal
# elements are those given in the first element of the pair.
# Thus, its columns are the eigenvectors corresponding to those eigenvalues.


Eigenfiddle:=proc(OurMat::Matrix,$)
    local i,n,real_eigens,eigenstuff,eigen_order;
  
  n:=RowDimension(OurMat);

  # The Maple function Eigenvectors returns a pair, the first of
  # which is a list of eigenvalues, and the second is a matrix
  # whose columns are the corresponding eigenvectors.
  # We ensure that the Matrix being processed is diagonal by
  # averaging it with its transpose.

  eigenstuff:=Eigenvectors(Matrix(n,n,(i,j)->(OurMat[i,j]+OurMat[j,i])/2,
                                    scan=diagonal[upper],shape=symmetric));

  # The following list contains pairs [eig,i], where i is the index
  # in the list. The idea is to sort the eigenvalues into increasing
  # order, but keep track of their original i's, so that we can
  # then use the same order for the eigenvectors.
 
  real_eigens:=[seq([eigenstuff[1][i],i],i=1..n)];

  # Now sort these into increasing values of the eigenvalues using
  # the pair_order function defined below.

  real_eigens:=sort(real_eigens,pair_order);

  # Get the index order - to be applied to the eigenvectors below.

  eigen_order:=map2(op,2,real_eigens);

  # Return pair,
  #   element 1 lists all (real) e-values
  #   element 2 is a transformation matrix, with the
  #              columns e-vectors of above e-values.

  [ map2(op,1,real_eigens), Matrix([Column(eigenstuff[2],eigen_order)]) ];
end:


# The following routine is used to order pairs in the routine Eigenfiddle().
# (the pairs are simply ordered according to their first elements).
# This procedure is only used within Eigenfiddle.

pair_order:=proc(eigenpair1::list(numeric),eigenpair2::list(numeric))
  evalb(eigenpair1[1]<eigenpair2[1]);
end:


# The following represents an arbitrary operator on the eigenspaces
# produced by the above procedure DigXspace_OpLC(): it uses three of
# the elements output by that routine - trying anything else is risky!
# The output is a "block matrix" of matrices, each of which represents
# the raw (adjusted SO(3)-reduced) transition amplitudes
#
#           <j_f,L_f || H || j_i,L_i>
#           -------------------------
#                sqrt (2*L_f+1)
#
# between the j_i state (Hamiltonian eigenstate) of AM L_i and the
# j_f state of AM L_f.
# Note that the (i1,i2) block matrix element corresponds to the
# AMs goodL[i1] and goodL[i2].
# This output is probably best displayed using the Show_Rats() procedure,
# (which is controlled by various global parameters), below.

# The argument tran_op will often be the quadrupole operator
# specified in quad_op, but anything else could be used here.
# The results need to be multiplied by sqrt(2L_f+1) to get
# the genuine SO(3)-reduced matrix elements. To get transition
# rate, we also need to multiply by 1/sqrt(2L_i+1). These factors
# are taken care of in the displaying function Show_Rats() below.

# (Note that the inverse of the eigenvector matrix (in each L space)
#  is required, and calculated within this function. This is probably
#  not inefficient because this function will usually only be
#  called once for a given set of parameters.)


AmpXspeig_OpLC:=proc(tran_op::list, eigen_bases::list,
                                           Xparams::list, goodL::list)
      local i,j,LL,L_min,L_max,L_factor,L_dims,L_ends,L_count,eigen_invs,
            tran_mat,block_tran_mat,ebase_tran_mat;

  L_count:=nops(goodL);

  if Lcount=0 then return fi:   # nothing to do

  L_min:=goodL[1]:
  L_max:=goodL[L_count]:

  # For the adjusted SO(3)-reduced transition operator, form the
  # transformation matrix encompassing all L values
  # (we cut it into blocks below).

  tran_mat:=RepXspace_OpLC(tran_op,op(Xparams),L_min,L_max);

  # Obtain the sizes of the blocks (one for each good L value).

  L_dims:=[seq(dimXspace(op(3..6,Xparams),LL),LL in goodL)];
  L_ends:=[seq(dimXspace(op(3..6,Xparams),L_min,LL),LL in goodL)];

  # Form the block matrix, each element of which is itself a matrix.
  # The (i,j) block is of size L_dims[i] x L_dims[j].

  block_tran_mat:=Matrix(L_count,
        (i,j)->SubMatrix(tran_mat,[L_ends[i]-L_dims[i]+1..L_ends[i]],
                                  [L_ends[j]-L_dims[j]+1..L_ends[j]]) );

  # Here we simply transform to the basis (an eigenbasis) specified
  # in the matrices eigen_basis. First form inverse transition matrices.

  eigen_invs:=map(MatrixInverse,eigen_bases);

  ebase_tran_mat:=Matrix(L_count,
        (i,j)->eigen_invs[i].block_tran_mat[i,j].eigen_bases[j]);

  # output block matrix

  ebase_tran_mat:

end;


###########################################################################
####------ Formatted output of eigenvalues and transition rates -------####
###########################################################################

# The following takes the 1st and 4th items of the quartet output from
# DigXspace_OpLC() and uses them to print out the eigenenergies in a
# nicely formatted list. Each value is taken with respect to the
# lowest eigenenergy, and scaled by the value of glb_eig_sft
# (which may be set using the procedure ACM_set_output() above,
# or through using ACM_Adapt() below).
# The value of toshow indicates the maximum number of eigenvalues
# to be shown at each value of the AM (which are listed in goodL).
# That set by ACM_set_listln is used as a default.
# The overall minimal value of the eigenenergy is returned.


Show_Eigs:=proc(eigen_vals::list,goodL::list,
                     toshow::nonnegint:=glb_eig_num,
                     L_min::nonnegint:=0, L_max::nonnegint:=1000000,$)
      local LL,i,eigen_low,L_count,L_top;
      global glb_low_pre,glb_rel_wid,glb_rel_pre,glb_eig_sft;

  if toshow=0 or nops(eigen_vals)=0 then return NULL fi:

  if nargs=4 then
    L_top:=L_min:
  else
    L_top:=L_max:
  fi:

  L_count:=nops(goodL):

  # Count how many L to output in range

  i:=0: 
  for LL to L_count do
    if goodL[LL]>=L_min and goodL[LL]<=L_top then
      i:=i+1:
    fi:
  od:

  if i=0 then
    return NULL
  fi:

  # Find smallest eigenvalue across all L spaces (hope all are real!)
  # Each sublist should be increasing.

  eigen_low:=min_head(eigen_vals);

  # display all required eigenvalues, with scaling given by glb_eig_sft.

  printf("Lowest eigenvalue is %.*f. Relative eigenvalues follow"
             " (each divided by %.*f):\n",
                       glb_low_pre,eigen_low,glb_low_pre,glb_eig_sft);

  for LL to L_count do
    if goodL[LL]>=L_min and goodL[LL]<=L_top then
      # print L and first relative eigenvalue
      printf("  At L=%2d: [%*.*f",goodL[LL],glb_rel_wid,glb_rel_pre,
                        (eigen_vals[LL][1]-eigen_low)/glb_eig_sft);
      # print remaining eigenvalues
      for i from 2 to min(nops(eigen_vals[LL]),toshow) do
        printf(",%*.*f",glb_rel_wid,glb_rel_pre,
                        (eigen_vals[LL][i]-eigen_low)/glb_eig_sft);
      od:
      # finish writing line
      printf("]\n");
    fi:
  od:

  eigen_low;   # output smallest eigenvalue - it might be needed!
end;


# The following routines are used to obtain, within a list of lists,
# the minimal value amongst the first elements.
# These procedures are only used within Show_Eigs().

min_head:=(alist)->min(op(map(fsel,alist)));
fsel:=(nlist)->`if`(nops(nlist)>0,nlist[1],NULL);


# The routine Show_Rats below prints out the selection of transition
# rates (and amplitudes) that are determined by showlist.
# Here goodL is the list of pertinent values of AM.
# These values label the blocks of the block matrix, the elements
# of which give the adjusted SO(3)-reduced matrix elements between
# the states having those AM values.

# Each element of showlist is either
#   1. A quartet of the form [Li,Lf,Li_index,Lf_index]
#      signifying a particular element to be output.
#      That value is the (Lf_index,Li_index) element of
#      the (Lf,Li) matrix in the block matrix block_mat.
#   2. A triple of the form [Li,Lf,Lf_index].
#      Then a vector of the above values is formed from
#      all possible values of Li_index (note strange order).
#   3. A pair of the form [Li,Lf].
#      Then a list of vectors of the above form is displayed
#      for all possible values of Lf_index.
#   4. A quintet of the form [Li,Lf,Li_index,Lf_index,L_mod].
#      This produces all quartets [LiX,LfX,Li_index,Lf_index]
#      with LiX=Li+k*L_mod and LfX=Lf+k*L_mod for k>=0.

# The argument toshow limits the number of values shown in
# the lists produced by [Li,Lf,Lf_index] and [Li,Lf].

# The argument amp_flag determines for which transition rates,
# the amplitudes are also displayed. If greater than 5, none are
# displayed. Otherwise, the amplitudes are displayed only for
# those format specifiers that have at least amp_flag elements.

# Each transition rate/amplitude is scaled by the value of
# glb_rat_sft or glb_amp_sft respectively
# (these may be set using the procedure ACM_set_output() above,
# or through using ACM_Adapt() below).

# In addition, the argument amp_mul_fun provides an additional
# scaling factor for the amplitudes, that depends on the
# AM values between which the operator acts. Thus amp_mul_fun
# is a procedure that takes two arguments Li and Lf.

# Note that the default values are determined each time the
# procedure is invoked, and not when it is initially defined.

# Also note that if the 3rd argument is not a list of lists,
# or the 6th or 7th arguments not procedure sthen
# an error will result (with that argument taken to be the 8th
# and an excess of arguments being flagged).

Show_Rats:=proc(block_mat::Matrix, goodL::list,
                  showlist::list(list):=glb_rat_lst,
                  toshow::integer:=glb_rat_num,
                  amp_flag::integer:=glb_amp_flg,
                  amp_fun::procedure:=quad_amp_fun,
                  rat_fun::procedure:=quad_rat_fun,$)
      local L1,L2,idx1,idx2,L1_off,L2_off,rate_ent,TR_matrix,
            TR_cols,TR_rows,Lmod,Lcount,col_count,
            item_preformat,item_format,
            rat_format1,rat_format2,amp_format1,amp_format2;
      global glb_low_pre,glb_rel_wid,glb_rel_pre,
             glb_rat_sft,glb_amp_sft,
             glb_rat_format1,glb_rat_format2,
             glb_amp_format1,glb_amp_format2:


  if nops(showlist)=0 then return fi:

  if ColumnDimension(block_mat)=0 then
    error "No transition rates avaiable!"
  fi:

  printf("Selected transition rates follow"
         " (each divided by %.*f", glb_low_pre,glb_rat_sft);

  if amp_flag<6 then
    printf(", amplitudes divided by %.*f", glb_low_pre,glb_amp_sft);
  fi:

  printf("):\n");

  # Specify format for the printing of each transition rate/amplitude.
  # Two stages - first sets the width and precision for values to output.

  item_preformat:= "%%%d.%df":
  item_format:=sprintf(item_preformat,glb_rel_wid,glb_rel_pre):

  # Change some %s specifications to 'integer' %d specifications.
  # (first two display single values, next two horizontal lists).

  rat_format1:=sprintf(glb_rat_format1,"%d","%d","%d","%d",item_format):
  amp_format1:=sprintf(glb_amp_format1,item_format):
  rat_format2:=sprintf(glb_rat_format2,"%d","%d","%d","%s"):
  amp_format2:=glb_amp_format2:

  # The scalings used for calculating the scaled transition rates
  # and the scaled amplitudes (which ought to be compatible) is
  # determined by the global variables glb_rat_sft and glb_amp_sft,
  # as used to define glb_rat_mul and glb_amp_mul below.
  # See section 6 of [Rowe2004], for useful info.
  # (Note that the dimSO3(L2) factor is required to convert to the
  # genuine SO(3)-reduced amplitude.)
  # The block_mat are multiplied by the glb_amp_mul to get the amplitude,
  # and after squaring, by glb_rat_mul to get the transition rate.
  # (The output is actually performed by the Show_Rat_element() and
  # Show_Amp_element() procedures below).

  # Now output the transition rates listed in showlist.
  # Note that those in the list that are not in the
  # range of those calculated are silently ignored.

  for rate_ent in showlist do

    if nops(rate_ent)<2 or nops(rate_ent)>5 then
      printf("  Bad transition rate specification: %a\n",rate_ent):
      next:
    fi:

    L1:=rate_ent[1]:
    L2:=rate_ent[2]:

    if L1<0 or L2<0 then next fi:

    if nops(rate_ent)<5 then

      # Locate indices in goodL for these Ls.

      if member(L1,goodL,'L1_off') and member(L2,goodL,'L2_off') then
        TR_matrix:=block_mat[L2_off,L1_off];
        TR_cols:=ColumnDimension(TR_matrix);
        TR_rows:=RowDimension(TR_matrix);

        if nops(rate_ent)=4 then    # output 4-index specifiers

          idx1:=rate_ent[3]:
          idx2:=rate_ent[4]:

          if idx1>0 and idx2>0 and idx1<=TR_cols and idx2<=TR_rows then
              printf(rat_format1,L1,idx1,L2,idx2,
                       evalf(rat_fun(L1,L2,TR_matrix[idx2,idx1])/glb_rat_sft)):

              if amp_flag<5 then
                printf(amp_format1,
                       evalf(amp_fun(L1,L2,TR_matrix[idx2,idx1])/glb_amp_sft)):
              fi:

              printf("\n"):
          fi:

        elif nops(rate_ent)=3 then    # output 3-index specifiers

          idx2:=rate_ent[3]:

          if idx2>0 and idx2<=TR_rows and TR_cols>0 then
              col_count:=min(TR_cols,toshow):

              printf(rat_format2,L1,L2,idx2,cat("[",
                  sprintf(item_format,
                      evalf(rat_fun(L1,L2,TR_matrix[idx2,1])/glb_rat_sft)),
                    seq(cat(",",
                      sprintf(item_format,
                         evalf(rat_fun(L1,L2,TR_matrix[idx2,i])/glb_rat_sft))),
                                        i=2..col_count),"]")):

              printf("\n"):

              if amp_flag<4 then
                printf(amp_format2,cat("[",
                  sprintf(item_format,
                      evalf(amp_fun(L1,L2,TR_matrix[idx2,1])/glb_amp_sft)),
                    seq(cat(",",
                      sprintf(item_format,
                         evalf(amp_fun(L1,L2,TR_matrix[idx2,i])/glb_amp_sft))),
                                        i=2..col_count),"]")):

                printf("\n"):
              fi:
          fi:

        else                          # output 2-index specifiers

          if TR_cols>0 and TR_rows>0 then
              col_count:=min(TR_cols,toshow):

              for idx2 to min(TR_rows,toshow) do
                printf(rat_format2,L1,L2,idx2,cat("[",
                  sprintf(item_format,
                      evalf(rat_fun(L1,L2,TR_matrix[idx2,1])/glb_rat_sft)),
                    seq(cat(",",
                      sprintf(item_format,
                         evalf(rat_fun(L1,L2,TR_matrix[idx2,i])/glb_rat_sft))),
                                        i=2..col_count),"]")):

                printf("\n"):

                if amp_flag<3 then
                printf(amp_format2,cat("[",
                  sprintf(item_format,
                      evalf(amp_fun(L1,L2,TR_matrix[idx2,1])/glb_amp_sft)),
                    seq(cat(",",
                      sprintf(item_format,
                         evalf(amp_fun(L1,L2,TR_matrix[idx2,i])/glb_amp_sft))),
                                        i=2..col_count),"]")):

                  printf("\n"):
                fi:
              od:
          fi:
        fi:
      fi:

    else       # output 5-index specifiers

        idx1:=rate_ent[3]:
        idx2:=rate_ent[4]:
        if idx1<=0 or idx2<=0 then next fi:
        Lmod:=rate_ent[5]:

        if Lmod>0 then    # Put Lcount+1 as number of rates required
          Lcount:=iquo(goodL[-1]-max(L1,L2),Lmod):
        elif Lmod<0 then
          Lmod:=-Lmod:
          Lcount:=iquo(min(L1,L2),Lmod):
          L1:=L1-Lcount*Lmod:
          L2:=L2-Lcount*Lmod:
        else
          Lcount:=0:
        fi:

        while Lcount>=0 do    # loop through all Lcount+1 cases

          if member(L1,goodL,'L1_off') and member(L2,goodL,'L2_off') then

            TR_matrix:=block_mat[L2_off,L1_off];
            TR_cols:=ColumnDimension(TR_matrix);
            TR_rows:=RowDimension(TR_matrix);

            if idx1<=TR_cols and idx2<=TR_rows then
                printf(rat_format1,L1,idx1,L2,idx2,
                       evalf(rat_fun(L1,L2,TR_matrix[idx2,idx1])/glb_rat_sft)):

                if amp_flag<6 then
                  printf(amp_format1,
                       evalf(amp_fun(L1,L2,TR_matrix[idx2,idx1])/glb_amp_sft)):
                fi:

                printf("\n"):
            fi:
          fi:

          L1:=L1+Lmod:
          L2:=L2+Lmod:
          Lcount:=Lcount-1;
        od

    fi:
  od:
  NULL;
end;

# A few routines to display individual transition rates and amplitudes,
# or lists of such. Each is multiplied by the value of glb_rat_mul or
# glb_amp_mul, which are defined (and vary with Li,Lf) in the above
# procedure.

#Get_Rat_element:=proc(TR_ME)
#  global glb_item_format,glb_rat_mul;
#  sprintf(glb_item_format,Re(TR_ME)^2*glb_rat_mul);
#end;
#
#Get_Amp_element:=proc(TR_ME)
#  global glb_item_format,glb_amp_mul;
#  sprintf(glb_item_format,Re(TR_ME)*glb_amp_mul);
#end;
#
#Get_Rat_col:=proc(TR_matrix,col_idx,row_count)
#  local i;
#
#  cat("[",Get_Rat_element(TR_matrix[1,col_idx]),
#      seq(cat(",",Get_Rat_element(TR_matrix[i,col_idx])),i=2..row_count),
#      "]"):
#end;
#
#Get_Amp_col:=proc(TR_matrix,col_idx,row_count)
#  local i;
#
#  cat("[",Get_Amp_element(TR_matrix[1,col_idx]),
#      seq(cat(",",Get_Amp_element(TR_matrix[i,col_idx])),i=2..row_count),
#      "]"):
#end;
#
#Get_Rat_row:=proc(TR_matrix,row_idx,col_count)
#  local i;
#
#  cat("[",Get_Rat_element(TR_matrix[row_idx,1]),
#      seq(cat(",",Get_Rat_element(TR_matrix[row_idx,i])),i=2..col_count),
#      "]"):
#end;
#
#Get_Amp_row:=proc(TR_matrix,row_idx,col_count)
#  local i;
#
#  cat("[",Get_Amp_element(TR_matrix[row_idx,1]),
#      seq(cat(",",Get_Amp_element(TR_matrix[row_idx,i])),i=2..col_count),
#      "]"):
#end;




# The routine Show_Mels below is an analogue of the above, but
# only prints out the adjusted matrix elements (amplitudes),
# Again, the values output are determined by showlist.
# Here goodL is the list of pertinent values of AM.
# These values label the blocks of the block matrix, the elements
# of which give the adjusted SO(3)-reduced matrix elements between
# the states having those AM values.

# Each element of showlist is either
#   1. A quartet of the form [Li,Lf,Li_index,Lf_index]
#      signifying a particular element to be output.
#      That value is the (Lf_index,Li_index) element of
#      the (Lf,Li) matrix in the block matrix block_mat.
#   2. A triple of the form [Li,Lf,Lf_index].
#      Then a vector of the above values is formed from
#      all possible values of Li_index (note strange order).
#   3. A pair of the form [Li,Lf].
#      Then a list of vectors of the above form is displayed
#      for all possible values of Lf_index.
#   4. A quintet of the form [Li,Lf,Li_index,Lf_index,L_mod].
#      This produces all quartets [LiX,LfX,Li_index,Lf_index]
#      with LiX=Li+k*L_mod and LfX=Lf+k*L_mod for k>=0.

# The argument toshow limits the number of values shown in
# the lists produced by [Li,Lf,Lf_index] and [Li,Lf].

# The argument amp_flag also serves to filter which elements
# from showlist get displayed (this is here for convenience -
# the same showlist can be used simultaneously here and in
# Show_Rats).
# If greater than 5, no matrix elements from showlist are displayed.
# Otherwise, matrix elements are displayed only for those
# format specifiers that have at least amp_flag elements.
# Default value is 0, meaning that all specifiers are enacted.

# Each transition rate/amplitude is scaled by the value of
# glb_rat_sft or glb_amp_sft respectively
# (these may be set using the procedure ACM_set_output() above,
# or through using ACM_Adapt() below).

# In addition, the argument amp_mul_fun provides an additional
# scaling factor for the matrix elements, that depends on the
# AM values between which the operator acts. Thus amp_mul_fun
# is a procedure that takes two arguments Li and Lf.
# (If omitted, the unit scaling is applied: i.e. it is 1).

# Note that the default values are determined each time the
# procedure is invoked, and not when it is initially defined.

# Also note that if the 3rd argument is not a list of lists,
# or the 6th argument is not a procedure then
# an error will result (with that argument taken to be the 7th
# and an excess of arguments being flagged).


Show_Mels:=proc(block_mat::Matrix, goodL::list,
                  showlist::list(list):=glb_rat_lst,
                  toshow::integer:=glb_rat_num,
                  amp_flag::integer:=0,
                  mel_fun::procedure:=mel_amp_fun,$)
      local L1,L2,idx1,idx2,L1_off,L2_off,rate_ent,TR_matrix,
            TR_cols,TR_rows,Lmod,Lcount,col_count,
            item_preformat,item_format,mel_format1,mel_format2;
      global glb_low_pre,glb_rel_wid,glb_rel_pre,
             glb_amp_sft,
             glb_mel_format1,glb_mel_format2:


  if nops(showlist)=0 then return fi:

  if ColumnDimension(block_mat)=0 then
    error "No matrix elements avaiable!"
  fi:

  printf("Selected matrix elements follow"
         " (each divided by %.*f):\n", glb_low_pre,glb_amp_sft);

  # Specify format for the printing of each transition rate/amplitude.
  # Two stages - first sets the width and precision for values to output.

  item_preformat:= "%%%d.%df":
  item_format:=sprintf(item_preformat,glb_rel_wid,glb_rel_pre):

  # Change some %s specifications to 'integer' %d specifications.

  mel_format1:=sprintf(glb_mel_format1,"%d","%d","%d","%d",item_format):
  mel_format2:=sprintf(glb_mel_format2,"%d","%d","%d","%s"):

  # The scalings used for calculating the scaled matrix elements
  # is determined by the global variable glb_amp_sft,
  # as used to define glb_amp_mul below.
  # The block_mat are multiplied by the glb_amp_mul to get the amplitude.
  # (The output is actually performed by the Show_Amp_element()
  # procedures below).

  # Now run through showlist.
  # Note that those in the list that are not in the range of those
  # present in (the pre-calculated) block_mat are silently ignored.

  for rate_ent in showlist do

    if nops(rate_ent)<2 or nops(rate_ent)>5 then
      printf("  Bad matrix element specification: %a\n",rate_ent):
      next:
    fi:

    L1:=rate_ent[1]:
    L2:=rate_ent[2]:

    if L1<0 or L2<0 then next fi:

    if nops(rate_ent)<5 then

      # Locate indices in goodL for these Ls.

      if member(L1,goodL,'L1_off') and member(L2,goodL,'L2_off') then
        TR_matrix:=block_mat[L2_off,L1_off];
        TR_cols:=ColumnDimension(TR_matrix);
        TR_rows:=RowDimension(TR_matrix);

        if nops(rate_ent)=4 and amp_flag<5 then   # output 4-index specifiers

          idx1:=rate_ent[3]:
          idx2:=rate_ent[4]:

          if idx1>0 and idx2>0 and idx1<=TR_cols and idx2<=TR_rows then

              printf(mel_format1,L1,idx1,L2,idx2,
                  evalf(mel_fun(L1,L2,TR_matrix[idx2,idx1])/glb_amp_sft)):

              printf("\n"):
          fi:

        elif nops(rate_ent)=3 and amp_flag<4 then # output 3-index specifiers

          idx2:=rate_ent[3]:

          if idx2>0 and idx2<=TR_rows and TR_cols>0 then
            col_count:=min(TR_cols,toshow):

            printf(mel_format2,L1,L2,idx2,cat("[",
                  sprintf(item_format,
                      evalf(mel_fun(L1,L2,TR_matrix[idx2,1])/glb_amp_sft)),
                    seq(cat(",",
                      sprintf(item_format,
                         evalf(mel_fun(L1,L2,TR_matrix[idx2,i])/glb_amp_sft))),
                                        i=2..col_count),"]")):
            printf("\n"):
          fi:

        elif nops(rate_ent)=2 and amp_flag<3 then # output 2-index specifiers

          if TR_cols>0 and TR_rows>0 then
            col_count:=min(TR_cols,toshow):

            for idx2 to min(TR_rows,toshow) do
              printf(mel_format2,L1,L2,idx2,cat("[",
                  sprintf(item_format,
                      evalf(mel_fun(L1,L2,TR_matrix[idx2,1])/glb_amp_sft)),
                    seq(cat(",",
                      sprintf(item_format,
                         evalf(mel_fun(L1,L2,TR_matrix[idx2,i])/glb_amp_sft))),
                                        i=2..col_count),"]")):

              printf("\n"):
            od:
          fi:
        fi:
      fi:

    elif amp_flag<6 then       # output 5-index specifiers

        idx1:=rate_ent[3]:
        idx2:=rate_ent[4]:
        if idx1<=0 or idx2<=0 then next fi:
        Lmod:=rate_ent[5]:

        if Lmod>0 then    # Put Lcount+1 as number of rates required
          Lcount:=iquo(goodL[-1]-max(L1,L2),Lmod):
        elif Lmod<0 then
          Lmod:=-Lmod:
          Lcount:=iquo(min(L1,L2),Lmod):
          L1:=L1-Lcount*Lmod:
          L2:=L2-Lcount*Lmod:
        else
          Lcount:=0:
        fi:

        while Lcount>=0 do    # loop through all Lcount+1 cases

          if member(L1,goodL,'L1_off') and member(L2,goodL,'L2_off') then

            TR_matrix:=block_mat[L2_off,L1_off];
            TR_cols:=ColumnDimension(TR_matrix);
            TR_rows:=RowDimension(TR_matrix);

            if idx1<=TR_cols and idx2<=TR_rows then
              printf(mel_format1,L1,idx1,L2,idx2,
                  evalf(mel_fun(L1,L2,TR_matrix[idx2,idx1])/glb_amp_sft)):
              printf("\n"):
            fi:
          fi:

          L1:=L1+Lmod:
          L2:=L2+Lmod:
          Lcount:=Lcount-1;
        od

    fi:
  od:
  NULL;
end;



###########################################################################

# The following generates a list of all pairs [L1,L2] in the
# range L_min .. L_max, and with | L1-L2 | <= L_diff, but with
# L1=1 and L2=1 omitted because such AM cannot occur).
# This can be used as input to ACM_add_rat_lst and ACM_set_rat_lst.

Lpairs_Gen:=proc(L_min::nonnegint,L_max::nonnegint,L_diff::nonnegint:=2,$)
  local pairL,L1,L2;

  if L_min>L_max then error "Nonsensical L range!" fi:

  if L_min=0 then
      pairL:=[[0,0],seq('[L,0],[0,L]',L=2..min(L_max,L_diff))]:
  else
    pairL:=[]:
  fi:

  for L2 from max(2,L_min) to L_max do
    pairL:=[op(pairL),seq([L1,L2],
                          L1=max(2,L_min,L2-L_diff)..min(L_max,L2+L_diff))]:
  od:

  pairL;
end;


###########################################################################
####--------------------- All purpose procedures ----------------------####
###########################################################################


# The following two procedures combine all the above in the easiest way.
# ACM_Scale() first diagonalises the Hamiltonian ham_op
# (using DigXspace_OpLC above)
# and displays the resulting eigenenergies per L value,
# relative to their overall lowest value (which is displayed),
# and the relative values are scaled by the global parameter glb_eig_sft.
# The output format is affected by the global parameters
# glb_eig_num, glb_rel_pre, glb_rel_wid, glb_low_pre.
# Then the amplitudes specified in the global list glb_amp_lst are
# printed, and then the transition rates specified in the global
# list glb_amp_lst are printed, both of the latter using the function
# Show_Rats.
# Everything displayed is scaled by global values as set above.

# The two procedures might never need their return values used. However,
# the return value is the triple
#                 [e_vals,Melements,goodL],
# where e_vals is a list of lists of eigenvalues (one list for each
# L-space in goodL), and Melements are the adjusted reduced matrix
# elements of transition rates (calculated in AmpXspeig_OpLC) stored
# as a block matrix. This latter is only calculated when the list of
#transition rates (glb_rat_lst) is non-empty (otherwise, tran_mat=NULL).
# The first and third argument may be used as the first two arguments
# to Show_Eigs() and the second and third as the first two argument to
# Show_Rats() to display further eigenenergies,
# transition rates/amplitudes without the need for recalculation.


ACM_Scale:=proc(ham_op::list,
                 anorm::algebraic, lambda_base::algebraic,
                 nu_min::nonnegint, nu_max::nonnegint,
                 v_min::nonnegint, v_max::nonnegint,
                 L_min::nonnegint, L_max::nonnegint,$)

  ACM_ScaleOrAdapt(0,0,args):
end;


# The following procedure performs as the above ACM_Scale(), 
# except that the global scaling parameters are determined within
# and then used.
# The parameters are set according to ACM_set_fits() routines etc.

# Note that this procedure sets the global scaling parameters.

ACM_Adapt:=proc(ham_op::list,
                 anorm::algebraic, lambda_base::algebraic,
                 nu_min::nonnegint,nu_max::nonnegint,
                 v_min::nonnegint,v_max::nonnegint,
                 L_min::nonnegint,L_max::nonnegint,$)

  ACM_ScaleOrAdapt(1,1,args):
end;


# The following procedure does the dirty work for the above two.

ACM_ScaleOrAdapt:=proc(fit_eig::nonnegint,fit_rat::nonnegint,
                     ham_op::list,
                     anorm::algebraic, lambda_base::algebraic,
                     nu_min::nonnegint, nu_max::nonnegint,
                     v_min::nonnegint, v_max::nonnegint,
                     L_min::nonnegint, L_max::nonnegint,$)
      local eigen_quin,tran_mat,goodL,eigen_low,L_mx,L1_off,L2_off;
      global glb_eig_num, quad_op, glb_rat_lst, glb_amp_lst, glb_amp_flg,
             glb_eig_fit, glb_eig_sft,
             glb_eig_L, glb_eig_idx,
             glb_rat_TRop, glb_amp_fun, glb_rat_fun,glb_amp_sft_fun,
             glb_rat_num, glb_rat_fit, glb_rat_sft, glb_amp_sft,
             glb_rat_L1, glb_rat_1dx, glb_rat_L2, glb_rat_2dx;

  if nargs=10 then  # no L_max
    L_mx:=L_min;
  else
    L_mx:=L_max
  fi:

  # When fitting values (if either fit_eig or fit_rat is non_zero)
  # we must ensure that the eigenvalue or transition rate with
  # respect to which we fit, and thus choose scaling parameters,
  # will actually be obtained. If not, we exit with an error message.
  # Note that we only perform this check when values of each variety
  # will actually be output (glb_eig_num>0 and nops(glb_rat_lst)>0 resp.)

  # First check the eigenenergy parameters

  if fit_eig>0 and glb_eig_num>0 then
    if glb_eig_L<L_min or glb_eig_L>L_mx or
         glb_eig_idx>dimXspace(nu_min,nu_max,v_min,v_max,glb_eig_L) then
      error "Reference state %1(%2) not available", glb_eig_L, glb_eig_idx;
    fi:
  fi:

  # Now check the parameters for the transition rates

  if fit_rat>0 and nops(glb_rat_lst)>0 then
    if glb_rat_L1<L_min or glb_rat_L1>L_mx or
         glb_rat_1dx>dimXspace(nu_min,nu_max,v_min,v_max,glb_rat_L1) then
      error "Reference state %1(%2) not available", glb_rat_L1, glb_rat_1dx;
    fi:

    if glb_rat_L2<L_min or glb_rat_L2>L_mx or
         glb_rat_2dx>dimXspace(nu_min,nu_max,v_min,v_max,glb_rat_L2) then
      error "Reference state %1(%2) not available", glb_rat_L2, glb_rat_2dx;
    fi:
  fi:

  # diagonalise the Hamiltonian on the specified space.
  # output is [ eigenval_list, goodL, Ps, Xparams ].

  eigen_quin:=DigXspace_OpLC(ham_op,args[4..-1]):
  goodL:=eigen_quin[4]:

  if glb_eig_num>0 then    # require eigenvalue output

    eigen_low:=min_head(eigen_quin[1]);   # determine smallest eigenvalue

    if fit_eig>0 then   # determine global scale factor for eigenvalues
      member(glb_eig_L,goodL,'LL'):   # find index for required L
      glb_eig_sft:=(eigen_quin[1][LL][glb_eig_idx]-eigen_low)/glb_eig_fit:

      if glb_eig_sft=0 then
        error "Cannot scale: reference state %1(%2) has lowest energy",
                                                    glb_eig_L, glb_eig_idx;
      fi
    fi:

    # display all required eigenvalues, with scaling given by glb_eig_sft.

    Show_Eigs(eigen_quin[1],goodL,glb_eig_num):
  fi:

  # Now turn attention to the transition rates, if any are required...

  if nops(glb_rat_lst)>0 then

    # obtain raw transition amplitudes

    tran_mat:=AmpXspeig_OpLC(glb_rat_TRop,op(2..-1,eigen_quin)):

    if fit_rat>0 then   # determine global scale factor for transition rates
      member(glb_rat_L1,goodL,'L1_off'):   # find indices for required Ls
      member(glb_rat_L2,goodL,'L2_off'):

      glb_rat_sft:=glb_rat_fun(glb_rat_L1,glb_rat_L2,
             tran_mat[L2_off,L1_off][glb_rat_2dx,glb_rat_1dx])/glb_rat_fit;

      if glb_rat_sft=0 then
        error "Cannot scale zero transition rate B(E2: %1(%2) -> %3(%4))",
                   glb_rat_L1,glb_rat_1dx,glb_rat_L2,glb_rat_2dx,glb_rat_fit;
      fi:

      # and set scaling factor for amplitudes

      glb_amp_sft:=glb_amp_sft_fun(glb_rat_sft);
    fi:

    # display required transition rates with scaling factor given
    # by glb_rat_sft, and that for amplitudes by glb_amp_sft.

    Show_Rats(tran_mat,goodL,glb_rat_lst,
                glb_rat_num,glb_amp_flg,glb_amp_fun,glb_rat_fun):

    # return the raw data in case more are required.

  else
    tran_mat:=Matrix(0,0):   # NULL matrix - to indicate that we have no MEs 
  fi:

  [eigen_quin[1],tran_mat,goodL]:

end;



###########################################################################
###########################################################################

# Some additional ACM stuff pertaining to the paper [RWC2009];
# the equation numbers are those in that paper.

# The following function creates Hamiltonians of the form given in (75).

# For the Hamiltonian of (89), set c1=1-2*alpha and c2=alpha
# (See also eqns. (4.230) & (4.220) of [RowanWood].)


RWC_ham:=(B,c1,c2,chi,kappa)->
  Hamiltonian_OpLC(-1/2/B,0,B*c1/2,B*c2/2,0,-chi,0,0,0,kappa);


# The following is its expectation value, given by (76).

RWC_expt:=proc(B::numeric,c1::numeric,c2::numeric,kappa::numeric,
                          a::numeric,lambda0::numeric,$)
  local aa:

  aa:=a^2:
  aa*(4+9/(lambda0-1))/8/B + B*lambda0*c1/2/aa
              + B*lambda0*(lambda0+1)*c2/2/aa^2 + kappa/2:
end:

# The following is its expectation value, given by (76), assuming that
# lambda is given in terms of aa by (82).

RWC_expt_link:=proc(B::numeric,c1::numeric,c2::numeric,kappa::numeric,
                          a::numeric,$)
  RWC_expt(args,evalf(RWC_dav(c1,c2,a))):
end:

RWC_dav:=proc(c1::numeric,c2::numeric,a::numeric,$)
  if c1>=0 then
    return 2.5:
  else
    return 1+sqrt(9+((a^2)*c1/c2)^2)/2:
  fi;
end:


# For the Davidson value of lambda given by
#    lambda(v)=1+sqrt( (v+3/2)^2 + C ),
# the following function returns a function of v that finds the
# nearest integer to
#    lambda(v)-lambda(0)
# which is of the same parity as v.
# This can be used as a basis type, for it is assured that
#    lambda(v+1)-lambda(v) = +1 or -1.

lambda_davi_fun:=proc(C::numeric)
  local difffun:

    difffun:=proc(v) option operator,arrow;
       local diffint:
       diffint:=floor(sqrt((v+1.5)^2+C)-sqrt(2.25+C)):
       if type(diffint-v,odd) then
         return diffint+1:
       else
         return diffint:
       fi:
    end:

   difffun:
end:

###########################################################################

# The procedure RWC_alam obtains values of the ACM parameters lambda
# and a, which are "optimal" in the cases of these Hamiltonians.
# This requires the minimal value of RWC_expt, given above.
# RWC1 and RWC2 are derivatives of RWC_expt for c1>0 and c1<0 resp
# (but scaled somewhat), with muf the value of 2(lambda-1) given by (82).
# [In deriving these results, I used d(mu)/d(aa)=aa*(c1/c2)^2/mu.]

# This procedure returns [a,lambda] as the optimal pair.

RWC_alam:=proc(B::numeric,c1::numeric,c2::numeric,$)
  local RWC1,RWC2,muf,aa0;

  if c1<0 then # (lambda is a function of aa)
               # there is always one positive solution here
               # (in fact, I've never found other real solns.)

    muf:=(aa) -> sqrt( 9+ (aa*c1/c2)^2 ):

    # The following is the derivative of the corrected (76) in [RWC2009].
    # But multiplied by 4*A^3*B*mu.
    RWC2:=(aa,mu) -> (c1/c2)^2 * (-9*aa^5/mu^2
                                    + aa^3*B^2*c1 + aa^2*B^2*c2*(mu+3))
                       + aa^3*(2*mu+9)
                       - B^2*mu*(mu+2)*(aa*c1+c2*(mu+4)):
    # The following results from the erroneous (76) in [RWC2009].
    #    RWC2:=(aa,mu) -> (c1/c2)^2 * (-9*aa^5/mu^2
    #                                  + 2*aa^3*B^2*c1 + 2*aa^2*B^2*c2*(mu+3))
    #                       + aa^3*(2*mu+9)
    #                       - 2*B^2*mu*(mu+2)*(aa*c1+c2*(mu+4)):
    aa0:=max(fsolve(RWC2(A,muf(A))=0,A)):
    return [sqrt(aa0),1+muf(aa0)/2]:


  else   # (lambda is constant here) fsolve produces real solns.
         # here there is always 1 positive soln;
         # and possibly two others that are negative.
         # Use max to exclude them.

    # corrected version
    RWC1:=(aa) -> aa^3 - B^2*c1*aa - 7*B^2*c2:
    # erroneous version
    # RWC1:=(aa) -> aa^3 - 2*B^2*c1*aa - 14*B^2*c2:

    aa0:=max(fsolve(RWC1(A)=0,A)):
    return [sqrt(aa0),2.5]:

  fi:

end:

# As above, but with lambda taken to be constant value of 2.5

RWC_alam_clam:=proc(B::numeric,c1::numeric,c2::numeric,$)
  local RWC1,RWC2,muf,aa0;

    RWC1:=(aa) -> aa^3 - B^2*c1*aa - 7*B^2*c2:
    aa0:=max(fsolve(RWC1(A)=0,A)):
    return [sqrt(aa0),2.5]:

end:


# This is DJR's original algorithm. It returns [a,lambda].

RWC_alam_old:=proc(B::numeric,alpha::numeric,step::numeric:=0.1)
  local lambda,ham_exp,x0,x1,H0,H1;

  if alpha<=0.5 then
    lambda:=proc(aa) option operator,arrow; 2.5 end
  else
    lambda:=proc(aa) option operator,arrow;
              1+sqrt(9+(aa*(2 - 1/alpha))^2)/2
            end
  fi;

  ham_exp:=proc(aa,lam) optionoperator,arrow;
         1/2*(aa*(lam+1.25))/(B*(lam - 1))
         +1/2*(B*(1 - 2*alpha)*lam)/aa
         +1/2*(B*alpha*lam*(lam+1))/aa^2
       end proc;

  if alpha<=0.2 then
    x0:=B
  elif alpha<=0.5 then
    x0:=B*sqrt(1.1 - 2*alpha)
  elif alpha>0.5 then
    x0:=B*sqrt(alpha - 0.5)
  fi:

  x1:=x0+step:
  H0:=ham_exp(x0,lambda(x0)); H1:=ham_exp(x1,lambda(x1));

  while H1<H0 do
    x0:=x1: H0:=H1:
    x1:=x1+step: H1:=ham_exp(x1,lambda(x1)):
  od;

  [sqrt(x0),lambda(x0)]
end:


# This does the same as RWC_alam above, but returns a third
# value, which is the name of a function of v that gives the (optimal)
# value of lambda(v)-lambda(0), an integer of same parity as v.


RWC_alam_fun:=proc(B::numeric,c1::numeric,c2::numeric,$)
  local RWC1,RWC2,muf,aa0;

  if c1<0 then # (lambda is a function of aa)
               # there is always one positive solution here
               # (in fact, I've never found other real solns.)

    muf:=(aa) -> sqrt( 9+ (aa*c1/c2)^2 ):

    # The following is the derivative of the corrected (76) in [RWC2009].
    # But multiplied by 4*A^3*B*mu.
    RWC2:=(aa,mu) -> (c1/c2)^2 * (-9*aa^5/mu^2
                                    + aa^3*B^2*c1 + aa^2*B^2*c2*(mu+3))
                       + aa^3*(2*mu+9)
                       - B^2*mu*(mu+2)*(aa*c1+c2*(mu+4)):
    aa0:=max(fsolve(RWC2(A,muf(A))=0,A)):
    return [sqrt(aa0),1+muf(aa0)/2,lambda_davi_fun((aa0*c1/c2/2)^2)]:


  else   # (lambda is constant here) fsolve produces real solns.
         # here there is always 1 positive soln;
         # and possibly two others that are negative.
         # Use max to exclude them.

    RWC1:=(aa) -> aa^3 - B^2*c1*aa - 7*B^2*c2:
    aa0:=max(fsolve(RWC1(A)=0,A)):
    return [sqrt(aa0),2.5,lambda_sho_fun]:

  fi:

end:



###########################################################################
# Finally, we read in the user-defined settings - including, most
# importantly, the location of the SO(5)>SO(3) Clebsch-Gordon coefficients.
###########################################################################

read "acm_user.mpl":

###########################################################################
# End of file.
###########################################################################

