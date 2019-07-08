# acm_user.mpl (version 1.4, 27 September 2015).

# This file contains user-definable settings for use by the
# Maple program acm1.4.mpl


# First, specify where to find the files containing the
# SO(5)>SO(3) CG coefficients (see below for a test).
# If desired, this value could be reset at the start of a worksheet.

SO5CG_directory:="/home/trevor/progs/so5/data/so5cg-data/":
#SO5CG_directory:="/home/twelsh1/projects/toronto/acm/so5cg-data/":
#SO5CG_directory:="/home/hs/staff/twelsh1/projects/toronto/acm/so5cg-data/":

# The SO5 CG filenames are of the form SO5CG_v1_v2-a2-L2_v3.
# Each of these is assumed to lie in the directory named
# SO5CG_v1_v2_v3 which itself lies in a directory
# named "v2=1/", "v2=2/", "v2=3/", ... , "v2=6/"
# below the directory specified above by SO5CG_directory.
# (Note that the final "/" in the directory name is necessary.)

# When uncommented, the following line tests the above directory,
# and (somewhat) the data therein (it should output two values,
# 0.522,0.431).

# show_CG_file(2,3,1,0,5):

# The following sets how lambda depends on seniority for the basis.
# 2 gives the ACM parity basis (0.0 is a placefiller, the third
# argument makes it silent).

ACM_set_basis_type(2,0.0,0):

# ACM_set_output() sets
#  the maximal number of decimal digits for output of values,
#  the maximal overall number of digits for output of values,
#  the maximal number of decimal digits for the lowest eigenvalue.
# (The fourth argument 0 makes it silent).

ACM_set_output(2,8,5,0):

# ACM_set_datum() specifies whether absolute eigenvalues are displayed (0),
# or are displayed relative to the minimal value (1).
# (The second argument 0 makes it silent).

ACM_set_datum(1,0):

# ACM_set_listln() sets
#  the maximal number of eigenvalues to show for each L,
#  the maximal number of transition rates in a "list",
# (The third argument 0 makes it silent).

ACM_set_listln(4,4,0):

# ACM_set_scales() sets the values of the factors which are used to
# divide the values of the (relative) eigenenergies and transition
# rates in the routine ACM_Scale().
# (The third argument 0 makes it silent).

ACM_set_scales(1.0,1.0,0):


# The procedures ACM_set_eig_fit() and ACM_set_rat_fit() specify
# values that are used by the procedure ACM_Adapt() to adjust
# the scaling factors so that a certain eigenvalue and a certain
# quadrupole transition rate take on specific values.
# Invoking ACM_set_eig_fit(val,L,ix) sets the eigenvalue
# scaling factor so that the ix'th lowest relative eigenvalue for
# angular momentum L takes the value val.
# (The fourth argument 0 makes it silent).

ACM_set_eig_fit(6.0,2,1,0):

# Invoking ACM_set_rat_fit(val,L1,L2,ix1,ix2) sets the
# quadrupole transition rate scalings so that the transition rate
# from the ix1'th lowest state of angular momentum L1 to the
# the ix2'th lowest state of angular momentum L2 takes the value val.
# (The sixth argument 0 makes it silent).

ACM_set_rat_fit(100.0,2,0,1,1,0):

# The transition rates that are to be displayed are stored in
# a list that is initialised using ACM_set_rat_lst(); this list
# can be augmented later using ACM_add_rat_lst().
# Here, we initially specify that no transition rates are to be displayed.

ACM_set_rat_lst( [] ):
ACM_set_amp_lst( [] ):

# The procedure ACM_amp_show(n) indicates for which elements of
# the list of transition rates, are the amplitudes to be displayed
# as well. 6 for none, otherwise only those designators with at
# least n parameters specified.

ACM_set_amp_show(2,0):

# Specify that transition rates and amplitudes will be calculated
# for the quadrupole operator.
# Specify that amplitude scale factors are sqrt of transition rate
# scale factors. And specify the format of output.

ACM_set_transition(quad_op,0):
ACM_set_sft_fun(sqrt_fun,0):

ACM_set_rat_form(quad_rat_fun,def_rat_format,def_rat_desg,0):
ACM_set_amp_form(quad_amp_fun,def_amp_format,def_amp_desg,0):

 
######################################################################

# Specify that when Maple chooses, it outputs 3 decimal places:

interface(displayprecision=3);

# Specify that matrices up to 15 x 15 get fully displyed:

interface(rtablesize=15);

# Use following for worksheet output:

#interface(prettyprint=3);

######################################################################
###  End of user-specified settings file acm1.4-user.mpl .
######################################################################

