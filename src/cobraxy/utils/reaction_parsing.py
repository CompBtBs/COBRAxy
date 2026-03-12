"""
Helpers to parse reaction strings into structured dictionaries.

Features:
- Reaction direction detection (forward, backward, reversible)
- Parsing of custom reaction strings into stoichiometric maps
- Conversion of a dict of raw reactions into a directional reactions dict
- Loading custom reactions from a tabular file (TSV)
"""
from enum import Enum
from typing import Dict
import re

try:
    from . import general_utils as utils
except:
    import general_utils as utils

# Reaction direction encoding:
class ReactionDir(Enum):
  """
  A reaction can go forward, backward, or be reversible (both directions).
  Cobrapy-style formulas encode direction using specific arrows handled here.
  """
  FORWARD    = "-->"
  BACKWARD   = "<--"
  REVERSIBLE = "<=>"

  @classmethod
  def fromReaction(cls, reaction :str) -> 'ReactionDir':
    """
    Takes a whole reaction formula string and looks for one of the arrows, returning the
    corresponding reaction direction.

    Args:
      reaction : the reaction's formula.
    
    Raises:
      ValueError : if no valid arrow is found.
    
    Returns:
      ReactionDir : the corresponding reaction direction.
    """
    for member in cls:
      if member.value in reaction: return member

    raise ValueError("No valid arrow found within reaction string.")

ReactionsDict = Dict[str, Dict[str, float]]


def add_custom_reaction(reactionsDict :ReactionsDict, rId :str, reaction :str) -> None:
  """
  Add one reaction entry to reactionsDict.

  The entry maps each substrate ID to its stoichiometric coefficient.
  If a substrate appears without an explicit coefficient, 1.0 is assumed.

  Args:
    reactionsDict: Dict to update in place.
    rId: Unique reaction ID.
    reaction: Reaction formula string.
  
  Returns:
    None

  Side effects: updates reactionsDict in place.
  """
  reaction = reaction.strip()
  if not reaction: return

  reactionsDict[rId] = {}
  # Assumes ' + ' is spaced to avoid confusion with charge symbols.
  for word in reaction.split(" + "):
    metabId, stoichCoeff = word, 1.0
    # Coefficient can be integer or float (dot decimal) and must be space-separated.
    foundCoeff = re.search(r"\d+(\.\d+)? ", word)
    if foundCoeff:
      wholeMatch  = foundCoeff.group(0)
      metabId     = word[len(wholeMatch):].strip()
      stoichCoeff = float(wholeMatch.strip())

    reactionsDict[rId][metabId] = stoichCoeff

  if not reactionsDict[rId]: del reactionsDict[rId] # Empty reactions are removed.


def create_reaction_dict(unparsed_reactions: Dict[str, str]) -> ReactionsDict:
    """
  Parse a dict of raw reaction strings into a directional reactions dict.

    Args:
    unparsed_reactions: Mapping reaction ID -> raw reaction string.

    Returns:
    ReactionsDict: Parsed dict. Reversible reactions produce two entries with _F and _B suffixes.
    """
    reactionsDict :ReactionsDict = {}
    for rId, reaction in unparsed_reactions.items():
        reactionDir = ReactionDir.fromReaction(reaction)
        left, right = reaction.split(f" {reactionDir.value} ")

    # Reversible reactions are split into two: forward (_F) and backward (_B).
        reactionIsReversible = reactionDir is ReactionDir.REVERSIBLE
        if reactionDir is not ReactionDir.BACKWARD:
            add_custom_reaction(reactionsDict, rId + "_F" * reactionIsReversible, left)
        
        if reactionDir is not ReactionDir.FORWARD:
            add_custom_reaction(reactionsDict, rId + "_B" * reactionIsReversible, right)
    
    return reactionsDict


def parse_custom_reactions(customReactionsPath :str) -> ReactionsDict:
  """
  Load custom reactions from a tabular file and parse into a reactions dict.

  Args:
    customReactionsPath: Path to the reactions file (TSV or CSV-like).
  
  Returns:
    ReactionsDict: Parsed reactions dictionary.
  """
  try:
    rows = utils.readCsv(utils.FilePath.fromStrPath(customReactionsPath), delimiter = "\t", skipHeader=False)
    if len(rows) <= 1:
      raise ValueError("The custom reactions file must contain at least one reaction.")

    id_idx, idx_formula = utils.findIdxByName(rows[0], "Formula")

  except Exception as e:
    # Fallback re-read with same settings; preserves original behavior
    rows = utils.readCsv(utils.FilePath.fromStrPath(customReactionsPath), delimiter = "\t", skipHeader=False)
    if len(rows) <= 1:
      raise ValueError("The custom reactions file must contain at least one reaction.")
    
    id_idx, idx_formula = utils.findIdxByName(rows[0], "Formula")
  
  reactionsData = {row[id_idx] : row[idx_formula] for row in rows[1:]}
  
  return create_reaction_dict(reactionsData)

