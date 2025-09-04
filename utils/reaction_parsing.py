from enum import Enum
import utils.general_utils as utils
from typing import Dict
import csv
import re

# Reaction direction encoding:
class ReactionDir(Enum):
  """
  A reaction can go forwards, backwards or be reversible (able to proceed in both directions).
  Models created / managed with cobrapy encode this information within the reaction's
  formula using the arrows this enum keeps as values.
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
  Adds an entry to the given reactionsDict. Each entry consists of a given unique reaction id
  (key) and a :dict (value) matching each substrate in the reaction to its stoichiometric coefficient.
  Keys and values are both obtained from the reaction's formula: if a substrate (custom metabolite id)
  appears without an explicit coeff, the value 1.0 will be used instead.

  Args:
    reactionsDict : dictionary encoding custom reactions information.
    rId : unique reaction id.
    reaction : the reaction's formula.
  
  Returns:
    None

  Side effects:
    reactionsDict : mut
  """
  reaction = reaction.strip()
  if not reaction: return

  reactionsDict[rId] = {}
  # We assume the '+' separating consecutive metabs in a reaction is spaced from them,
  # to avoid confusing it for electrical charge:
  for word in reaction.split(" + "):
    metabId, stoichCoeff = word, 1.0
    # Implicit stoichiometric coeff is equal to 1, some coeffs are floats.

    # Accepted coeffs can be integer or floats with a dot (.) decimal separator
    # and must be separated from the metab with a space:
    foundCoeff = re.search(r"\d+(\.\d+)? ", word)
    if foundCoeff:
      wholeMatch  = foundCoeff.group(0)
      metabId     = word[len(wholeMatch):].strip()
      stoichCoeff = float(wholeMatch.strip())

    reactionsDict[rId][metabId] = stoichCoeff

  if not reactionsDict[rId]: del reactionsDict[rId] # Empty reactions are removed.


def create_reaction_dict(unparsed_reactions: Dict[str, str]) -> ReactionsDict:
    """
    Parses the given dictionary into the correct format.

    Args:
        unparsed_reactions (Dict[str, str]): A dictionary where keys are reaction IDs and values are unparsed reaction strings.

    Returns:
        ReactionsDict: The correctly parsed dict.
    """
    reactionsDict :ReactionsDict = {}
    for rId, reaction in unparsed_reactions.items():
        reactionDir = ReactionDir.fromReaction(reaction)
        left, right = reaction.split(f" {reactionDir.value} ")

        # Reversible reactions are split into distinct reactions, one for each direction.
        # In general we only care about substrates, the product information is lost.
        reactionIsReversible = reactionDir is ReactionDir.REVERSIBLE
        if reactionDir is not ReactionDir.BACKWARD:
            add_custom_reaction(reactionsDict, rId + "_F" * reactionIsReversible, left)
        
        if reactionDir is not ReactionDir.FORWARD:
            add_custom_reaction(reactionsDict, rId + "_B" * reactionIsReversible, right)
        
        # ^^^ to further clarify: if a reaction is NOT reversible it will not be marked as _F or _B
        # and whichever direction we DO keep (forward if --> and backward if <--) loses this information.
        # This IS a small problem when coloring the map in marea.py because the arrow IDs in the map follow
        # through with a similar convention on ALL reactions and correctly encode direction based on their
        # model of origin. TODO: a proposed solution is to unify the standard in RPS to fully mimic the maps,
        # which involves re-writing the "reactions" dictionary.
    
    return reactionsDict


def parse_custom_reactions(customReactionsPath :str) -> ReactionsDict:
  """
  Creates a custom dictionary encoding reactions information from a csv file containing
  data about these reactions, the path of which is given as input.

  Args:
    customReactionsPath : path to the reactions information file.
  
  Returns:
    ReactionsDict : dictionary encoding custom reactions information.
  """
  reactionsData :Dict[str, str] = {row[0]: row[1] for row in utils.readCsv(utils.FilePath.fromStrPath(customReactionsPath), delimiter = "\t")} 
  return create_reaction_dict(reactionsData)

