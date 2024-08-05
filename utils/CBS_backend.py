from swiglpk import *
import random
import pandas as pd
import numpy as np
import cobra as cb

# Initialize LP problem
def initialize_lp_problem(S):

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
    
    

# Solve LP problem from the structure of the metabolic model
def create_and_solve_lp_problem(lb,ub,nrows, ncol, len_vector, ia, ja, ar, 
                                obj_coefs,reactions,return_lp=False):
    
    
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
    
    
# Solve LP problem from the structure of the metabolic model
def solve_lp_problem(lp,obj_coefs,reactions):
   
    # Set the coefficients of the objective function
    i=1
    for ind_coef in obj_coefs:
        glp_set_obj_coef(lp, i, ind_coef);
        i+=1

    # Initialize the parameters    
    params=glp_smcp()
    params.presolve=GLP_ON
    params.msg_lev = GLP_MSG_ALL
    params.tm_lim=4000
    glp_init_smcp(params)
    
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
    

# Create LP structure
def create_lp_structure(model):
    
    reactions=[el.id for el in model.reactions]
    coefs_obj=[reaction.objective_coefficient for reaction in model.reactions]
    
    # Lower and upper bounds
    lb=[reaction.lower_bound for reaction in model.reactions]
    ub=[reaction.upper_bound for reaction in model.reactions]
    
    # Create S matrix
    S=cb.util.create_stoichiometric_matrix(model,array_type="dok")
    
    return S,lb,ub,coefs_obj,reactions

# CBS sampling interface
def randomObjectiveFunctionSampling(model, nsample, coefficients_df, df_sample):

    S,lb,ub,coefs_obj,reactions = create_lp_structure(model)
    len_vector, values, indexes, ia, ja, ar, nrow, ncol = initialize_lp_problem(S)
    
    for i in range(nsample):
      
        coefs_obj=coefficients_df.iloc[:,i].values
            
        if coefs_obj[-1]==1: #minimize
            coefs_obj= coefs_obj[0:-1] * -1
        else:
            coefs_obj=coefs_obj[0:-1]

        fluxes,Z = create_and_solve_lp_problem(lb,ub, nrow, ncol, len_vector, 
                                                        ia, ja, ar, coefs_obj,reactions,return_lp=False)
        df_sample.loc[i] = fluxes 
    pass

def randomObjectiveFunctionSampling_cobrapy(model, nsample, coefficients_df, df_sample):
    
    for i in range(nsample):

        dict_coeff={}
        if(coefficients_df.iloc[-1][i]==1):
            type_problem = -1 #minimize
        else:
            type_problem = 1
            
        for rxn in [reaction.id for reaction in model.reactions]:
            dict_coeff[model.reactions.get_by_id(rxn)] = coefficients_df.loc[rxn][i] * type_problem
            
        model.objective = dict_coeff
        solution =  model.optimize().fluxes
        for rxn, flux in solution.items():
            df_sample.loc[i][rxn] = flux

    pass

# Create random coefficients for CBS
def randomObjectiveFunction(model, n_samples, df_fva, seed=0):
    

        #reactions = model.reactions
        reactions = [reaction.id for reaction in model.reactions]
        cont=seed
        list_ex=reactions.copy()
        list_ex.append("type_of_problem")
        coefficients_df = pd.DataFrame(index=list_ex,columns=[str(i) for i in range(n_samples)])

        for i in range(0, n_samples):
           
            cont=cont+1
            random.seed(cont)
            
            # Genera un numero casuale tra 0 e 1
            threshold = random.random() #coefficiente tra 0 e 1
            
            for reaction in reactions:

                cont=cont+1
                random.seed(cont)
                        
                val=random.random()   
                
                if val>threshold:

                    cont=cont+1
                    random.seed(cont)                           
                   
                    c=2*random.random()-1 #coefficiente tra -1 e 1
                    
                    val_max=np.max([df_fva.loc[reaction,"minimum"],df_fva.loc[reaction,"maximum"]])
                    
                    if val_max!=0: #solo se la fva è diversa da zero
                        coefficients_df.loc[reaction,str(i)] = c/val_max #divido per la fva
                    else:
                        coefficients_df.loc[reaction,str(i)] = 0

                else:
                    coefficients_df.loc[reaction,str(i)] = 0

            cont=cont+1
            random.seed(cont)
                    
            if random.random()<0.5:
                coefficients_df.loc["type_of_problem",str(i)] = 0 #maximize
            else:
                coefficients_df.loc["type_of_problem",str(i)] = 1 #minimize
                            
        return coefficients_df