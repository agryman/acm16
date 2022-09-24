# generate data files from Example 4.5 for Python tests
read "acm1.4a-examples-startup.mpl":

# Make the directory if it doesn't exist.
make_dir := proc(dir)
    if not FileTools[Exists](dir) then
        mkdir(dir):
    end if
end:

# Generate RepXspace() test data.
generate_repxspace := proc(base_dir, B, RWC_ham)
    local dir, rep:

    dir := String(base_dir, "repxspace/"):
    make_dir(dir):

    rep := RepXspace(RWC_ham, sqrt(B), 2.5, 0, 5, 0, 5, 0):

    ExportMatrix(String(dir, "rwc-ham.csv"), rep):
end:

# Generate RepXspace_Prod() test data.
generate_repxspace_prod := proc(base_dir, B)
    local dir, cases, case, rep, path:

    dir := String(base_dir, "repxspace-prod/"):
    make_dir(dir):

    cases:=[
        [[Radial_D2b],"radial-d2b"],
        [[Radial_bm2],"radial-bm2"],
        [[Radial_b],"radial-b"],
        [[Radial_b2],"radial-b2"],
        [[SpHarm_310],"spharm-310"]
    ]:

    for case in cases do
        rep := RepXspace_Prod(case[1], sqrt(B), 2.5, 0, 5, 0, 5, 0, 0):
        path := String(dir, case[2], ".csv"):
        ExportMatrix(path, rep):
    end do:
end:

# Generate RepSO5r3_Prod_rem() test data.
generate_repso5r3_prod_rem := proc(base_dir)
    local dir, rep:

    dir := String(base_dir, "repso5r3-prod-rem/"):
    make_dir(dir):

    rep := RepSO5r3_Prod_rem([SpHarm_310], 0, 5, 0, 0):

    ExportMatrix(String(dir, "spharm-310.csv"), rep):
end:

# Generate all test data for the specified basis type.
generate := proc(basis_type)
    local base_dir, B, c1, c2, chi, kappa, x1, x3, x4, x6, x10, RWC_ham:

    ACM_set_basis_type(basis_type);
    ACM_show_lambda_fun();

    make_dir("generated/"):
    make_dir("generated/example-4-5/"):
    make_dir("generated/example-4-5/05050/");

    base_dir := String("generated/example-4-5/05050/", basis_type, "/"):
    make_dir(base_dir):

    # Example 4.1 Hamiltonian
    B:=50: c2:=2.0: c1:=1-2*c2: chi:=1.5: kappa:=1.0:
    x1:=-1/2/B: x3:=B*c1/2: x4:=B*c2/2: x6:=-chi: x10:=kappa:
    RWC_ham:=ACM_Hamiltonian(x1,0,x3,x4,0,x6,0,0,0,x10);

    generate_repxspace(base_dir, B, RWC_ham):
    generate_repxspace_prod(base_dir, B):
    generate_repso5r3_prod_rem(base_dir):
end:
