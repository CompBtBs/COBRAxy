"""
CBS backend utilities using GLPK for constraint-based sampling.

This module builds and solves LP problems from COBRA models, supports random
objective function generation (CBS), and provides both swiglpk-based and
COBRApy-based sampling fallbacks.
"""
from swiglpk import *
import random
import pandas as pd
import numpy as np
import cobra as cb

# Initialize LP problem
def initialize_lp_problem(S):
    """
    Prepare sparse matrix structures for GLPK given a stoichiometric matrix.

    Args:
        S: Stoichiometric matrix (DOK or similar) as returned by cobra.util.create_stoichiometric_matrix.

    Returns:
        tuple: (len_vector, values, indexes, ia, ja, ar, nrows, ncol)
            - len_vector: number of non-zero entries
            - values: list of non-zero values
            - indexes: list of (row, col) indices for non-zero entries
            - ia, ja, ar: GLPK-ready arrays for the sparse matrix
            - nrows, ncol: matrix shape
    """

    len_vector=len(S.keys())
    values=list(S.values())
    indexes=list(S.keys())
    ia = intArray(len_vector+1); 
    ja = intArray(len_vector+1);
    ar = doubleArray(len_vector+1);
    
    i=0
    ind_row=[indexes[i][0]+1 for i in range(0, len(values) )]
    ind_col=[indexes[i][1]+1 for i in range(0, len(values) )]
    for i in range(1, len(values) + 1): 
        ia[i]=ind_row[i-1]
        ja[i]=ind_col[i-1]
        ar[i] = values[i-1]
    
    nrows=S.shape[0]
    ncol=S.shape[1]
    
    return len_vector, values, indexes, ia, ja, ar, nrows, ncol
    
    

def create_and_solve_lp_problem(lb,ub,nrows, ncol, len_vector, ia, ja, ar, 
                                obj_coefs,reactions,return_lp=False):
    """
    Create and solve a GLPK LP problem for a metabolic model.

    Args:
        lb, ub: Lower/upper bounds per reaction (lists of floats).
        nrows, ncol: Dimensions of the S matrix.
        len_vector, ia, ja, ar: Sparse matrix data prepared for GLPK.
        obj_coefs: Objective function coefficients (list of floats).
        reactions: Reaction identifiers (list of str), used for output mapping.
        return_lp: If True, also return the GLPK problem object (caller must delete).

    Returns:
        tuple: (fluxes, Z) or (fluxes, Z, lp) if return_lp=True.
    """
    
    
    lp = glp_create_prob();
    glp_set_prob_name(lp, "sample");
    glp_set_obj_dir(lp, GLP_MAX);
    glp_add_rows(lp, nrows);
    eps = 1e-16
    for i in range(nrows):
        glp_set_row_name(lp, i+1, "constrain_"+str(i+1));
        glp_set_row_bnds(lp, i+1, GLP_FX, 0.0, 0.0);
    glp_add_cols(lp, ncol);
    for i in range(ncol):
        glp_set_col_name(lp, i+1, "flux_"+str(i+1));
        glp_set_col_bnds(lp, i+1, GLP_DB,lb[i]-eps,ub[i]+eps);
    glp_load_matrix(lp, len_vector, ia, ja, ar);
    
    try:
        fluxes,Z=solve_lp_problem(lp,obj_coefs,reactions)
        if return_lp:
            return fluxes,Z,lp
        else:
            glp_delete_prob(lp);
            return fluxes,Z
    except Exception as e:
        glp_delete_prob(lp)
        raise Exception(e)
    
    
def solve_lp_problem(lp,obj_coefs,reactions):
    """
    Configure and solve an LP with GLPK using provided objective coefficients.

    Args:
        lp: GLPK problem handle.
        obj_coefs: Objective coefficients.
        reactions: Reaction identifiers (unused in computation; length used for output).

    Returns:
        tuple: (fluxes, Z) where fluxes is a list of primal column values and Z is the objective value.
    """
   
    # Set the coefficients of the objective function
    i=1
    for ind_coef in obj_coefs:
        glp_set_obj_coef(lp, i, ind_coef);
        i+=1

    # Initialize the parameters    
    params=glp_smcp()
    params.presolve=GLP_ON
    params.msg_lev = GLP_MSG_ERR
    params.tm_lim=4000
    glp_init_smcp(params)

    glp_term_out(GLP_OFF)

    try:
    
        # Solve the problem
        glp_scale_prob(lp,GLP_SF_AUTO)
        
        value=glp_simplex(lp, params) 

        Z = glp_get_obj_val(lp);

        if value == 0:
            fluxes = []
            for i in range(len(reactions)): fluxes.append(glp_get_col_prim(lp, i+1))
            return fluxes,Z
        else:
            raise Exception("error in LP problem. Problem:",str(value)) 
    except Exception as e:
        # Re-enable terminal output for error reporting
        glp_term_out(GLP_ON)
        raise Exception(e)
    finally:
        # Re-enable terminal output after solving
        glp_term_out(GLP_ON)
    
def create_lp_structure(model):
    """
    Extract the LP structure (S matrix, bounds, objective) from a COBRA model.

    Args:
        model (cobra.Model): The COBRA model.

    Returns:
        tuple: (S, lb, ub, coefs_obj, reactions)
    """
    
    reactions=[el.id for el in model.reactions]
    coefs_obj=[reaction.objective_coefficient for reaction in model.reactions]
    
    # Lower and upper bounds
    lb=[reaction.lower_bound for reaction in model.reactions]
    ub=[reaction.upper_bound for reaction in model.reactions]
    
    # Create S matrix
    S=cb.util.create_stoichiometric_matrix(model,array_type="dok")
    
    return S,lb,ub,coefs_obj,reactions

def randomObjectiveFunctionSampling(model, nsample, coefficients_df, df_sample):
    """
    Sample fluxes using GLPK by iterating over random objective functions.

    Args:
        model: COBRA model.
        nsample: Number of samples to generate.
        coefficients_df: DataFrame of objective coefficients (columns are samples; last row is type_of_problem).
        df_sample: Output DataFrame to fill with flux samples (index: sample, columns: reactions).

    Returns:
        None
    """

    S,lb,ub,coefs_obj,reactions = create_lp_structure(model)
    len_vector, values, indexes, ia, ja, ar, nrow, ncol = initialize_lp_problem(S)
    
    for i in range(nsample):
      
        coefs_obj=coefficients_df.iloc[:,i].values
            
        if coefs_obj[-1]==1: # minimize
            coefs_obj= coefs_obj[0:-1] * -1
        else:
            coefs_obj=coefs_obj[0:-1]

        fluxes,Z = create_and_solve_lp_problem(lb,ub, nrow, ncol, len_vector, 
                                                        ia, ja, ar, coefs_obj,reactions,return_lp=False)
        df_sample.loc[i] = fluxes 
    return

def randomObjectiveFunctionSampling_cobrapy(model, nsample, coefficients_df, df_sample):
    """
    Fallback sampling using COBRApy's optimize with per-sample randomized objectives.

    Args:
        model: COBRA model.
        nsample: Number of samples to generate.
        coefficients_df: DataFrame of objective coefficients (columns are samples; last row is type_of_problem).
        df_sample: Output DataFrame to fill with flux samples (index: sample, columns: reactions).

    Returns:
        None
    """
    
    for i in range(nsample):

        dict_coeff={}
        if(coefficients_df.iloc[-1][i]==1):
            type_problem = -1 # minimize
        else:
            type_problem = 1 # maximize
            
        for rxn in [reaction.id for reaction in model.reactions]:
            dict_coeff[model.reactions.get_by_id(rxn)] = coefficients_df.loc[rxn][i] * type_problem
            
        model.objective = dict_coeff
        solution =  model.optimize().fluxes
        for rxn, flux in solution.items():
            df_sample.loc[i][rxn] = flux

    return

def randomObjectiveFunction(model, n_samples, df_fva, seed=0):
    """
    Create random objective function coefficients for CBS sampling.

    The last row 'type_of_problem' encodes 0 for maximize and 1 for minimize.

    Args:
        model: COBRA model.
        n_samples: Number of random objectives to generate.
        df_fva: Flux Variability Analysis results with 'minimum' and 'maximum' per reaction.
        seed: Seed for reproducibility.

    Returns:
        pandas.DataFrame: Coefficients DataFrame indexed by reaction IDs plus 'type_of_problem'.
    """
    # reactions = model.reactions
    reactions = [reaction.id for reaction in model.reactions]
    cont = seed
    list_ex = reactions.copy()
    list_ex.append("type_of_problem")
    coefficients_df = pd.DataFrame(index=list_ex, columns=[str(i) for i in range(n_samples)])

    for i in range(0, n_samples):

        cont = cont + 1
        random.seed(cont)

        # Generate a random threshold in [0, 1]
        threshold = random.random()

        for reaction in reactions:

            cont = cont + 1
            random.seed(cont)

            val = random.random()

            if val > threshold:

                cont = cont + 1
                random.seed(cont)

                # Coefficient in [-1, 1]
                c = 2 * random.random() - 1

                val_max = np.max([abs(df_fva.loc[reaction, "minimum"]), abs(df_fva.loc[reaction, "maximum"])])

                if val_max != 0: # only if FVA is non-zero
                    coefficients_df.loc[reaction, str(i)] = c / val_max # scale by FVA
                else:
                    coefficients_df.loc[reaction, str(i)] = 0

            else:
                coefficients_df.loc[reaction, str(i)] = 0

        cont = cont + 1
        random.seed(cont)

        if random.random() < 0.5:
            coefficients_df.loc["type_of_problem", str(i)] = 0 # maximize
        else:
            coefficients_df.loc["type_of_problem", str(i)] = 1 # minimize

    return coefficients_df
