import numpy as np
import numpy.linalg as linalg

def placeholder_norm_sq_operator(RHS_, placeholder_norm_sq_parameters):
    RHS_ = np.array(RHS_).reshape(-1, 1)
    tmp_weight_ = placeholder_norm_sq_parameters
    if tmp_weight_ is None:
        tmp_weight_ = np.ones_like(RHS_)
    else:
        tmp_weight_ = np.array(tmp_weight_).reshape(-1, 1)
    tmp_norm_sq = np.dot(RHS_.T, tmp_weight_*RHS_)
    return tmp_norm_sq

def placeholder_forward_operator(LHS_, placeholder_forward_parameters):
    LHS_ = np.array(LHS_).reshape(-1, 1)
    A__ = np.array(placeholder_forward_parameters)
    A_LHS_ = np.dot(A__, LHS_)
    return A_LHS_

def placeholder_adjoint_operator(RHS_, placeholder_adjoint_parameters):
    RHS_ = np.array(RHS_).reshape(-1, 1)
    A__ = np.array(placeholder_adjoint_parameters)
    AT_RHS_ = np.dot(A__.T, RHS_)
    return AT_RHS_

def placeholder_cg_lsq(
        flag_verbose,
        placeholder_norm_sq_operator,
        placeholder_norm_sq_parameters, 
        placeholder_forward_operator,
        placeholder_forward_parameters, 
        placeholder_adjoint_operator,
        placeholder_adjoint_parameters,
        img_stack_, 
        tolerance_cg_lsq =1e-5, 
        n_iteration=None
        ):
    str_thisfunction = 'placeholder_cg_lsq'
    if flag_verbose is None:
        flag_verbose = 0
    # RHS_ stores the right hand side of the linear system
    RHS_ = placeholder_adjoint_operator(
        img_stack_,
        placeholder_adjoint_parameters
        )
    if n_iteration is None:
        n_iteration = len(RHS_)
    if flag_verbose:
        print('entering: ' + str_thisfunction + ': n_iteration = ' + str(n_iteration) + ' tolerance_cg_lsq = ' + str(tolerance_cg_lsq))
    # LHS_ stores the current solution
    LHS_ = np.zeros_like(RHS_);
    # RES_ stores the current residual
    RES_ = RHS_ - placeholder_adjoint_operator(
            placeholder_forward_operator(
                LHS_,
                placeholder_forward_parameters
                ),
            placeholder_adjoint_parameters
            )
    # PCG_ stores the current search direction
    PCG_ = RES_ 
    beta_num = placeholder_norm_sq_operator(RES_,placeholder_norm_sq_parameters)
    beta_den = 1.0
    beta = beta_num / max(1e-12,beta_den)
    niteration = 0
    flag_continue = 1
    while flag_continue:
        if flag_verbose:
            print(' % niteration = ' + str(niteration) + ' beta_num = ' + str(beta_num))
        zeta = placeholder_norm_sq_operator(
            placeholder_forward_operator(
                PCG_,
                placeholder_forward_parameters
                ),
            placeholder_norm_sq_parameters)
        alph = beta_num / max(1e-12,zeta)
        if flag_verbose:
            print(' % beta_num = ' + str(beta_num))
            print(' % zeta = ' + str(zeta))
            print(' % alph = ' + str(alph))
        LHS_ = LHS_ + alph * PCG_
        RES_ = RES_ - alph * placeholder_adjoint_operator(
            placeholder_forward_operator(
                PCG_,
                placeholder_forward_parameters
                ),
            placeholder_adjoint_parameters)
        beta_den = beta_num
        beta_num = placeholder_norm_sq_operator(
            RES_,
            placeholder_norm_sq_parameters
            )
        beta = beta_num / max(1e-12,beta_den)
        PCG_ = RES_ + beta * PCG_
        niteration += 1
        if niteration >= n_iteration:
            flag_continue = 0
        if np.sqrt(beta_num) < tolerance_cg_lsq:
            flag_continue = 0
    if flag_verbose:
        print(str_thisfunction + ': niteration = ' + str(niteration) + ' beta_num = ' + str(beta_num))
    if flag_verbose:
        print('exiting: ' + str_thisfunction)
    return (LHS_,np.sqrt(beta_num),niteration)

def test_cg_lsq():
    flag_verbose = 1
    A__ = np.array([[1, 2 ,-1], [3, 4, -2], [5, 6, 3]])
    AtA__ = np.dot(A__.T, A__)
    if flag_verbose:
        print('AtA__ = ' + str(AtA__))  
    placeholder_forward_parameters = A__
    placeholder_adjoint_parameters = A__
    placeholder_norm_sq_parameters = None
    LHS_tru_ = np.array([[1], [2] , [3]])
    img_stack_ = placeholder_forward_operator(LHS_tru_, placeholder_forward_parameters)
    RHS_tru_ = placeholder_adjoint_operator(img_stack_, placeholder_adjoint_parameters)
    LHS_nrm_ = np.linalg.solve(AtA__,RHS_tru_)
    LHS_rec, tmp_err_, tmp_nitera_ = placeholder_cg_lsq(
        flag_verbose,
        placeholder_norm_sq_operator,
        placeholder_norm_sq_parameters,
        placeholder_forward_operator,
        placeholder_forward_parameters,
        placeholder_adjoint_operator,
        placeholder_adjoint_parameters,
        img_stack_
        )   
    if flag_verbose:
        print(' % LHS_tru_ vs LHS_rec_: ' + str(linalg.norm(LHS_tru_ - LHS_rec)/linalg.norm(LHS_tru_)))
 
if __name__ == '__main__':
    test_cg_lsq()
    print('Final error should be small as long as the normal matrix is well conditioned.')
