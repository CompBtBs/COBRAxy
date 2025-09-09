import os
import csv
import cobra
import pickle
import argparse
import pandas as pd
from typing import Optional, Tuple, Union, List, Dict
import utils.general_utils as utils
import utils.rule_parsing  as rulesUtils

################################- DATA GENERATION -################################
ReactionId = str
def generate_rules(model: cobra.Model, *, asParsed = True) -> Union[Dict[ReactionId, rulesUtils.OpList], Dict[ReactionId, str]]:
    """
    Generates a dictionary mapping reaction ids to rules from the model.

    Args:
        model : the model to derive data from.
        asParsed : if True parses the rules to an optimized runtime format, otherwise leaves them as strings.

    Returns:
        Dict[ReactionId, rulesUtils.OpList] : the generated dictionary of parsed rules.
        Dict[ReactionId, str] : the generated dictionary of raw rules.
    """
    # Is the below approach convoluted? yes
    # Ok but is it inefficient? probably
    # Ok but at least I don't have to repeat the check at every rule (I'm clinically insane)
    _ruleGetter   =  lambda reaction : reaction.gene_reaction_rule
    ruleExtractor = (lambda reaction :
        rulesUtils.parseRuleToNestedList(_ruleGetter(reaction))) if asParsed else _ruleGetter

    return {
        reaction.id : ruleExtractor(reaction)
        for reaction in model.reactions
        if reaction.gene_reaction_rule }

def generate_reactions(model :cobra.Model, *, asParsed = True) -> Dict[ReactionId, str]:
    """
    Generates a dictionary mapping reaction ids to reaction formulas from the model.

    Args:
        model : the model to derive data from.
        asParsed : if True parses the reactions to an optimized runtime format, otherwise leaves them as they are.

    Returns:
        Dict[ReactionId, str] : the generated dictionary.
    """

    unparsedReactions = {
        reaction.id : reaction.reaction
        for reaction in model.reactions
        if reaction.reaction 
    }

    if not asParsed: return unparsedReactions
    
    return reactionUtils.create_reaction_dict(unparsedReactions)

def get_medium(model:cobra.Model) -> pd.DataFrame:
    trueMedium=[]
    for r in model.reactions:
        positiveCoeff=0
        for m in r.metabolites:
            if r.get_coefficient(m.id)>0:
                positiveCoeff=1;
        if (positiveCoeff==0 and r.lower_bound<0):
            trueMedium.append(r.id)

    df_medium = pd.DataFrame()
    df_medium["reaction"] = trueMedium
    return df_medium

def generate_bounds(model:cobra.Model) -> pd.DataFrame:

    rxns = []
    for reaction in model.reactions:
        rxns.append(reaction.id)

    bounds = pd.DataFrame(columns = ["lower_bound", "upper_bound"], index=rxns)

    for reaction in model.reactions:
        bounds.loc[reaction.id] = [reaction.lower_bound, reaction.upper_bound]
    return bounds



def generate_compartments(model: cobra.Model) -> pd.DataFrame:
    """
    Generates a DataFrame containing compartment information for each reaction.
    Creates columns for each compartment position (Compartment_1, Compartment_2, etc.)
    
    Args:
        model: the COBRA model to extract compartment data from.
        
    Returns:
        pd.DataFrame: DataFrame with ReactionID and compartment columns
    """
    pathway_data = []

    # First pass: determine the maximum number of pathways any reaction has
    max_pathways = 0
    reaction_pathways = {}

    for reaction in model.reactions:
        # Get unique pathways from all metabolites in the reaction
        if type(reaction.annotation['pathways']) == list:
            reaction_pathways[reaction.id] = reaction.annotation['pathways']
            max_pathways = max(max_pathways, len(reaction.annotation['pathways']))
        else:
            reaction_pathways[reaction.id] = [reaction.annotation['pathways']]

    # Create column names for pathways
    pathway_columns = [f"Pathway_{i+1}" for i in range(max_pathways)]

    # Second pass: create the data
    for reaction_id, pathways in reaction_pathways.items():
        row = {"ReactionID": reaction_id}
        
        # Fill pathway columns
        for i in range(max_pathways):
            col_name = pathway_columns[i]
            if i < len(pathways):
                row[col_name] = pathways[i]
            else:
                row[col_name] = None  # or "" if you prefer empty strings

        pathway_data.append(row)

    return pd.DataFrame(pathway_data)