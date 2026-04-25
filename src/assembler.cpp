PetscErrorCode Assembler::assemble(Mat& K, Vec& f)
{
    const int ndof = static_cast<int>(mesh.num_nodes()) * 2;

    // Create PETSc matrix and vector
    PetscErrorCode ierr;
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, ndof, ndof, 12, nullptr, &K); CHKERRQ(ierr);
    ierr = MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE); CHKERRQ(ierr);
    ierr = MatSetFromOptions(K); CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, ndof, &f); CHKERRQ(ierr);
    ierr = VecZeroEntries(f); CHKERRQ(ierr);

    // Assemble element contributions
    for (size_t e = 0; e < mesh.num_elements(); ++e) {
        std::array<double, 36> ke;
        double area;
        element_stiffness_CST(static_cast<int>(e), ke, area);

        const auto& elem = mesh.element(e);

        // DOF indices: [2*v0, 2*v0+1, 2*v1, 2*v1+1, 2*v2, 2*v2+1]
        PetscInt dofs[6] = {
            2 * elem.vid[0],     2 * elem.vid[0] + 1,
            2 * elem.vid[1],     2 * elem.vid[1] + 1,
            2 * elem.vid[2],     2 * elem.vid[2] + 1
        };

        PetscScalar vals[36];
        for (int i = 0; i < 36; ++i) vals[i] = ke[i];

        ierr = MatSetValues(K, 6, dofs, 6, dofs, vals, ADD_VALUES); CHKERRQ(ierr);
    }

    ierr = MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    // Apply Dirichlet BCs via penalty method
    // Penalty value: large relative to max diagonal
    PetscScalar max_diag = 0.0;
    for (int i = 0; i < ndof; ++i) {
        PetscScalar val;
        MatGetValue(K, i, i, &val);
        if (val > max_diag) max_diag = val;
    }
    const PetscScalar penalty = max_diag * 1e14;

    for (const auto& bc : dbc) {
        const PetscInt row = static_cast<PetscInt>(bc.dof);
        // Zero row and column, put penalty on diagonal
        ierr = MatZeroRowsColumns(K, 1, &row, penalty, nullptr, nullptr); CHKERRQ(ierr);
        // Set RHS = penalty * prescribed_value
        ierr = VecSetValue(f, row, penalty * bc.val, INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(f); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(f); CHKERRQ(ierr);

    return 0;
}
