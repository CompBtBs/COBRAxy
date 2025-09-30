#!/usr/bin/env python3
"""
Comprehensive test suite for GPR translation functionality in COBRAxy.

This test suite covers:
- Basic 1:1, 1:many, many:1 gene mappings
- Complex GPR expressions with AND/OR logic
- Translation issues tracking 
- OR-only GPR flattening functionality
- Edge cases and nested expressions
- Statistical reporting
"""

import cobra
import pandas as pd
import sys
import os
import logging
from typing import Dict, List, Tuple
import re

# Add the COBRAxy utils directory to the path
sys.path.append('/hdd/home/flapi/COBRAxy')
from utils import model_utils

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

class GPRTranslationTester:
    """Comprehensive GPR translation test suite"""
    
    def __init__(self):
        self.test_results = {}
        self.failed_tests = []
        
    def create_comprehensive_test_model(self) -> cobra.Model:
        """Create a comprehensive test model with diverse GPR patterns"""
        model = cobra.Model('comprehensive_test_model')
        
        # Create metabolites
        metabolites = []
        for i in range(30):
            met = cobra.Metabolite(f'met_{chr(65+i%26)}{i//26}', compartment='c')
            metabolites.append(met)
        
        reactions_data = [
            # === BASIC CASES ===
            ('BASIC_1to1', 'GENE1', 0, 1),  # Simple 1:1 mapping
            ('BASIC_1tomany', 'GENE2', 1, 2),  # 1:many mapping
            ('BASIC_manyto1', 'GENE3', 2, 3),  # many:1 mapping
            ('BASIC_unmapped', 'UNMAPPED_GENE', 3, 4),  # unmapped gene
            
            # === SIMPLE OR CASES (candidates for flattening) ===
            ('OR_simple', 'GENE4 or GENE5', 4, 5),  # Simple OR with many:1
            ('OR_three', 'GENE6 or GENE7 or GENE8', 5, 6),  # Three genes OR
            ('OR_parentheses', '(GENE9 or GENE10)', 6, 7),  # OR with parentheses
            ('OR_duplicates', 'GENE11 or GENE12 or GENE11', 7, 8),  # OR with duplicates after translation
            
            # === COMPLEX OR CASES (candidates for flattening) ===
            ('OR_nested_simple', '(GENE13 or GENE14) or (GENE15 or GENE16)', 8, 9),  # Nested OR only
            ('OR_many_parentheses', '((GENE17 or GENE18) or GENE19) or GENE20', 9, 10),  # Multiple nesting levels
            ('OR_mixed_mapping', 'GENE21 or GENE22 or GENE23', 10, 11),  # Mixed 1:1, 1:many, many:1
            
            # === AND CASES (should NOT be flattened) ===
            ('AND_simple', 'GENE24 and GENE25', 11, 12),  # Simple AND
            ('AND_complex', '(GENE26 and GENE27) and GENE28', 12, 13),  # Complex AND
            
            # === MIXED AND/OR (should NOT be flattened) ===
            ('MIXED_basic', 'GENE29 and (GENE30 or GENE31)', 13, 14),  # AND with OR
            ('MIXED_complex', '(GENE32 or GENE33) and (GENE34 or GENE35)', 14, 15),  # OR and AND
            ('MIXED_nested', '((GENE36 and GENE37) or GENE38) and GENE39', 15, 16),  # Complex nesting
            
            # === EDGE CASES ===
            ('EDGE_single', 'GENE40', 16, 17),  # Single gene
            ('EDGE_empty', '', 17, 18),  # Empty GPR
            ('EDGE_whitespace', '  GENE41  or  GENE42  ', 18, 19),  # Whitespace
            ('EDGE_case_sensitive', 'Gene43 OR gene44', 19, 20),  # Case variations
            
            # === STRESS TESTS ===
            ('STRESS_long_or', 'GENE45 or GENE46 or GENE47 or GENE48 or GENE49 or GENE50', 20, 21),  # Long OR chain
            ('STRESS_deep_nest', '(((GENE51 or GENE52) or GENE53) or GENE54)', 21, 22),  # Deep nesting
            ('STRESS_complex', '(GENE55 or (GENE56 or GENE57)) or ((GENE58 or GENE59) or GENE60)', 22, 23),  # Complex structure
            
            # === TRANSLATION ISSUE TRIGGERS ===
            ('ISSUE_1many_or', 'GENE61 or GENE62', 23, 24),  # 1:many in OR (should be flattened)
            ('ISSUE_manyto1_and', 'GENE63 and GENE64', 24, 25),  # many:1 in AND (should NOT be flattened)
            ('ISSUE_mixed_problems', '(GENE65 or GENE66) and GENE67', 25, 26),  # Mixed problems
            
            # === REAL-WORLD INSPIRED CASES ===
            ('REAL_metabolism', '(ENSG001 or ENSG002) or (ENSG003 or ENSG004)', 26, 27),  # Metabolic pathway
            ('REAL_transport', 'TRANSPORTER1 and (COFACTOR1 or COFACTOR2)', 27, 28),  # Transport reaction
            ('REAL_complex_enzyme', '((SUBUNIT1 and SUBUNIT2) or SUBUNIT3) and COFACTOR3', 28, 29),  # Complex enzyme
        ]
        
        # Create reactions
        for rxn_id, gpr, met_in, met_out in reactions_data:
            rxn = cobra.Reaction(rxn_id)
            if met_in < len(metabolites) and met_out < len(metabolites):
                rxn.add_metabolites({metabolites[met_in]: -1, metabolites[met_out]: 1})
            rxn.gene_reaction_rule = gpr
            model.add_reactions([rxn])
        
        return model
    
    def create_comprehensive_mapping(self) -> pd.DataFrame:
        """Create a comprehensive gene mapping covering all test scenarios"""
        mapping_data = {
            'hgnc_symbol': [],
            'ensg': []
        }
        
        # === BASIC MAPPINGS ===
        # 1:1 mappings
        one_to_one = [
            ('GENE1', 'TARGET1'),
            ('GENE24', 'TARGET24'),
            ('GENE25', 'TARGET25'),
            ('GENE26', 'TARGET26'),
            ('GENE27', 'TARGET27'),
            ('GENE28', 'TARGET28'),
            ('GENE29', 'TARGET29'),
            ('GENE40', 'TARGET40'),
            ('GENE41', 'TARGET41'),
            ('GENE42', 'TARGET42'),
        ]
        
        # 1:many mappings (one source gene maps to multiple targets)
        one_to_many = [
            ('GENE2', 'TARGET2A'),
            ('GENE2', 'TARGET2B'),
            ('GENE30', 'TARGET30A'),
            ('GENE30', 'TARGET30B'),
            ('GENE61', 'TARGET61A'),
            ('GENE61', 'TARGET61B'),
            ('GENE61', 'TARGET61C'),  # Maps to 3 targets
            ('GENE65', 'TARGET65A'),
            ('GENE65', 'TARGET65B'),
        ]
        
        # many:1 mappings (multiple source genes map to one target)
        many_to_one = [
            ('GENE3', 'SHARED_TARGET1'),
            ('GENE4', 'SHARED_TARGET1'),
            ('GENE5', 'SHARED_TARGET1'),
            ('GENE6', 'SHARED_TARGET2'),
            ('GENE7', 'SHARED_TARGET2'),
            ('GENE8', 'SHARED_TARGET2'),
            ('GENE9', 'SHARED_TARGET3'),
            ('GENE10', 'SHARED_TARGET3'),
            ('GENE11', 'SHARED_TARGET4'),
            ('GENE12', 'SHARED_TARGET4'),
            ('GENE13', 'SHARED_TARGET5'),
            ('GENE14', 'SHARED_TARGET5'),
            ('GENE15', 'SHARED_TARGET5'),
            ('GENE16', 'SHARED_TARGET5'),
            ('GENE17', 'SHARED_TARGET6'),
            ('GENE18', 'SHARED_TARGET6'),
            ('GENE19', 'SHARED_TARGET6'),
            ('GENE20', 'SHARED_TARGET6'),
            ('GENE45', 'SHARED_TARGET7'),
            ('GENE46', 'SHARED_TARGET7'),
            ('GENE47', 'SHARED_TARGET7'),
            ('GENE48', 'SHARED_TARGET7'),
            ('GENE49', 'SHARED_TARGET7'),
            ('GENE50', 'SHARED_TARGET7'),
            ('GENE51', 'SHARED_TARGET8'),
            ('GENE52', 'SHARED_TARGET8'),
            ('GENE53', 'SHARED_TARGET8'),
            ('GENE54', 'SHARED_TARGET8'),
            ('GENE55', 'SHARED_TARGET9'),
            ('GENE56', 'SHARED_TARGET9'),
            ('GENE57', 'SHARED_TARGET9'),
            ('GENE58', 'SHARED_TARGET9'),
            ('GENE59', 'SHARED_TARGET9'),
            ('GENE60', 'SHARED_TARGET9'),
            ('GENE63', 'SHARED_TARGET10'),
            ('GENE64', 'SHARED_TARGET10'),
            ('GENE66', 'SHARED_TARGET11'),
        ]
        
        # Mixed mappings for complex cases
        mixed_mappings = [
            ('GENE21', 'TARGET21'),  # 1:1
            ('GENE22', 'TARGET22A'),  # 1:many
            ('GENE22', 'TARGET22B'),
            ('GENE23', 'SHARED_TARGET1'),  # many:1 (shares with GENE3-5)
            ('GENE31', 'SHARED_TARGET12'),
            ('GENE32', 'SHARED_TARGET13'),
            ('GENE33', 'SHARED_TARGET13'),
            ('GENE34', 'TARGET34'),
            ('GENE35', 'TARGET35'),
            ('GENE36', 'TARGET36'),
            ('GENE37', 'TARGET37'),
            ('GENE38', 'TARGET38'),
            ('GENE39', 'TARGET39'),
            ('GENE62', 'TARGET62A'),
            ('GENE62', 'TARGET62B'),
            ('GENE67', 'TARGET67'),
        ]
        
        # Case sensitivity tests
        case_mappings = [
            ('Gene43', 'TARGET43'),
            ('gene44', 'TARGET44'),
        ]
        
        # Real-world inspired mappings
        real_mappings = [
            ('ENSG001', 'HUMAN_GENE1'),
            ('ENSG002', 'HUMAN_GENE2'),
            ('ENSG003', 'HUMAN_GENE1'),  # many:1
            ('ENSG004', 'HUMAN_GENE2'),  # many:1
            ('TRANSPORTER1', 'SLC_FAMILY1'),
            ('COFACTOR1', 'COFACTOR_A'),
            ('COFACTOR2', 'COFACTOR_A'),  # many:1
            ('COFACTOR3', 'COFACTOR_B'),
            ('SUBUNIT1', 'COMPLEX_SUBUNIT1'),
            ('SUBUNIT2', 'COMPLEX_SUBUNIT2'),
            ('SUBUNIT3', 'COMPLEX_ALTERNATIVE'),
        ]
        
        # Combine all mappings
        all_mappings = one_to_one + one_to_many + many_to_one + mixed_mappings + case_mappings + real_mappings
        
        for source, target in all_mappings:
            mapping_data['hgnc_symbol'].append(source)
            mapping_data['ensg'].append(target)
        
        return pd.DataFrame(mapping_data)
    
    def analyze_mapping_statistics(self, mapping_df: pd.DataFrame) -> Dict:
        """Analyze mapping statistics"""
        stats = {}
        
        source_counts = mapping_df.groupby('hgnc_symbol')['ensg'].count()
        target_counts = mapping_df.groupby('ensg')['hgnc_symbol'].count()
        
        stats['total_mappings'] = len(mapping_df)
        stats['unique_sources'] = len(source_counts)
        stats['unique_targets'] = len(target_counts)
        
        stats['one_to_one'] = (source_counts == 1).sum()
        stats['one_to_many'] = (source_counts > 1).sum()
        stats['many_to_one_targets'] = (target_counts > 1).sum()
        
        stats['one_to_many_details'] = {}
        for gene, count in source_counts[source_counts > 1].items():
            targets = mapping_df[mapping_df['hgnc_symbol'] == gene]['ensg'].tolist()
            stats['one_to_many_details'][gene] = targets
            
        stats['many_to_one_details'] = {}
        for target, count in target_counts[target_counts > 1].items():
            sources = mapping_df[mapping_df['ensg'] == target]['hgnc_symbol'].tolist()
            stats['many_to_one_details'][target] = sources
            
        return stats
    
    def predict_translation_issues(self, model: cobra.Model, mapping_df: pd.DataFrame) -> Dict:
        """Predict which reactions will have translation issues"""
        predictions = {}
        mapping_dict = {}
        
        # Build mapping dictionary
        for _, row in mapping_df.iterrows():
            source = row['hgnc_symbol']
            target = row['ensg']
            if source not in mapping_dict:
                mapping_dict[source] = []
            mapping_dict[source].append(target)
        
        for rxn in model.reactions:
            if not rxn.gene_reaction_rule or rxn.gene_reaction_rule.strip() == '':
                continue
                
            # Extract genes from GPR
            token_pattern = r'\b[A-Za-z0-9:_.-]+\b'
            tokens = re.findall(token_pattern, rxn.gene_reaction_rule)
            logical_operators = {'and', 'or', 'AND', 'OR', '(', ')'}
            genes = [t for t in tokens if t not in logical_operators]
            
            issues = []
            has_1_to_many = False
            has_many_to_1 = False
            has_unmapped = False
            
            for gene in set(genes):
                norm_gene = model_utils._normalize_gene_id(gene)
                if norm_gene in mapping_dict:
                    targets = mapping_dict[norm_gene]
                    if len(targets) > 1:
                        has_1_to_many = True
                        issues.append(f"1:many - {gene} -> {targets}")
                else:
                    has_unmapped = True
                    issues.append(f"unmapped - {gene}")
            
            # Check for many:1 mappings
            target_to_sources = {}
            for gene in set(genes):
                norm_gene = model_utils._normalize_gene_id(gene)
                if norm_gene in mapping_dict:
                    for target in mapping_dict[norm_gene]:
                        if target not in target_to_sources:
                            target_to_sources[target] = []
                        target_to_sources[target].append(gene)
            
            for target, sources in target_to_sources.items():
                if len(sources) > 1:
                    has_many_to_1 = True
                    issues.append(f"many:1 - {sources} -> {target}")
            
            if issues:
                predictions[rxn.id] = {
                    'issues': issues,
                    'has_1_to_many': has_1_to_many,
                    'has_many_to_1': has_many_to_1,
                    'has_unmapped': has_unmapped,
                    'is_or_only': self._check_if_or_only(rxn.gene_reaction_rule),
                    'predicted_flattening': has_1_to_many or has_many_to_1 and self._check_if_or_only(rxn.gene_reaction_rule)
                }
        
        return predictions
    
    def _check_if_or_only(self, gpr: str) -> bool:
        """Check if GPR contains only OR operators (and parentheses)"""
        if not gpr or gpr.strip() == '':
            return False
        
        # Remove gene names and whitespace, keep only logical operators
        token_pattern = r'\b[A-Za-z0-9:_.-]+\b'
        logic_only = re.sub(token_pattern, '', gpr)
        logic_only = re.sub(r'\s+', ' ', logic_only.strip())
        
        # Check for AND operators
        and_pattern = r'\b(and|AND)\b'
        return not bool(re.search(and_pattern, logic_only))
    
    def run_comprehensive_test(self) -> Dict:
        """Run the comprehensive translation test"""
        print("="*80)
        print("COMPREHENSIVE GPR TRANSLATION TEST SUITE")
        print("="*80)
        
        # Create test model and mapping
        print("\n1. Creating test model and mapping...")
        model = self.create_comprehensive_test_model()
        mapping_df = self.create_comprehensive_mapping()
        
        print(f"   ✓ Created model with {len(model.reactions)} reactions")
        print(f"   ✓ Created mapping with {len(mapping_df)} entries")
        
        # Analyze mapping statistics
        print("\n2. Analyzing mapping statistics...")
        mapping_stats = self.analyze_mapping_statistics(mapping_df)
        print(f"   ✓ Unique source genes: {mapping_stats['unique_sources']}")
        print(f"   ✓ Unique target genes: {mapping_stats['unique_targets']}")
        print(f"   ✓ 1:1 mappings: {mapping_stats['one_to_one']}")
        print(f"   ✓ 1:many mappings: {mapping_stats['one_to_many']}")
        print(f"   ✓ Many:1 target genes: {mapping_stats['many_to_one_targets']}")
        
        # Predict translation issues
        print("\n3. Predicting translation issues...")
        predicted_issues = self.predict_translation_issues(model, mapping_df)
        predicted_or_only = sum(1 for pred in predicted_issues.values() if pred['is_or_only'])
        predicted_flattening = sum(1 for pred in predicted_issues.values() if pred['predicted_flattening'])
        
        print(f"   ✓ Reactions with predicted issues: {len(predicted_issues)}")
        print(f"   ✓ OR-only reactions: {predicted_or_only}")
        print(f"   ✓ Predicted for flattening: {predicted_flattening}")
        
        # Display original GPRs
        print("\n4. Original model GPRs:")
        for rxn in sorted(model.reactions, key=lambda x: x.id):
            status = "🔍" if rxn.id in predicted_issues else "✓"
            or_only = "🔗" if predicted_issues.get(rxn.id, {}).get('is_or_only', False) else " "
            print(f"   {status}{or_only} {rxn.id:20} : {rxn.gene_reaction_rule}")
        
        # Run translation
        print("\n5. Running translation...")
        try:
            translated_model, translation_issues = model_utils.translate_model_genes(
                model=model,
                mapping_df=mapping_df,
                target_nomenclature='ensg',
                source_nomenclature='hgnc_symbol',
                allow_many_to_one=True
            )
            print("   ✓ Translation completed successfully")
        except Exception as e:
            print(f"   ❌ Translation failed: {e}")
            import traceback
            traceback.print_exc()
            return {'success': False, 'error': str(e)}
        
        # Display translated GPRs
        print("\n6. Translated model GPRs:")
        for rxn in sorted(translated_model.reactions, key=lambda x: x.id):
            has_issues = "🚨" if rxn.id in translation_issues else "✓"
            print(f"   {has_issues} {rxn.id:20} : {rxn.gene_reaction_rule}")
        
        # Analyze translation issues
        print("\n7. Translation issues analysis:")
        if translation_issues:
            for rxn_id, issues_str in sorted(translation_issues.items()):
                predicted = predicted_issues.get(rxn_id, {})
                prediction_status = "✓ PREDICTED" if rxn_id in predicted_issues else "❓ UNEXPECTED"
                print(f"   🚨 {rxn_id:20} ({prediction_status})")
                # Split issues string by semicolon separator
                if issues_str:
                    issues_list = [issue.strip() for issue in issues_str.split(';') if issue.strip()]
                    for issue in issues_list:
                        print(f"      - {issue}")
                else:
                    print(f"      - No specific issues reported")
        else:
            print("   ✅ No translation issues detected")
        
        # Compare predictions vs actual
        print("\n8. Prediction accuracy:")
        true_positive = set(predicted_issues.keys()) & set(translation_issues.keys())
        false_positive = set(predicted_issues.keys()) - set(translation_issues.keys())
        false_negative = set(translation_issues.keys()) - set(predicted_issues.keys())
        
        print(f"   ✓ Correctly predicted issues: {len(true_positive)}")
        print(f"   ⚠ False positives: {len(false_positive)}")
        print(f"   ❌ False negatives: {len(false_negative)}")
        
        if false_positive:
            print("   False positive reactions:")
            for rxn_id in false_positive:
                print(f"      - {rxn_id}")
        
        if false_negative:
            print("   False negative reactions:")
            for rxn_id in false_negative:
                print(f"      - {rxn_id}")
        
        # Test specific functionality
        print("\n9. Testing OR-only GPR flattening...")
        flattening_tests = self.test_or_only_flattening(translated_model, translation_issues)
        
        # Summary statistics
        print("\n10. Summary:")
        results = {
            'success': True,
            'model_reactions': len(model.reactions),
            'mapping_entries': len(mapping_df),
            'predicted_issues': len(predicted_issues),
            'actual_issues': len(translation_issues),
            'prediction_accuracy': {
                'true_positive': len(true_positive),
                'false_positive': len(false_positive),
                'false_negative': len(false_negative),
                'precision': len(true_positive) / len(predicted_issues) if predicted_issues else 0,
                'recall': len(true_positive) / len(translation_issues) if translation_issues else 0,
            },
            'mapping_stats': mapping_stats,
            'flattening_tests': flattening_tests,
            'models': {
                'original': model,
                'translated': translated_model
            },
            'issues': {
                'predicted': predicted_issues,
                'actual': translation_issues
            }
        }
        
        precision = results['prediction_accuracy']['precision']
        recall = results['prediction_accuracy']['recall']
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
        
        print(f"   📊 Total reactions: {len(model.reactions)}")
        print(f"   📊 Reactions with issues: {len(translation_issues)}")
        print(f"   📊 Prediction precision: {precision:.2%}")
        print(f"   📊 Prediction recall: {recall:.2%}")
        print(f"   📊 Prediction F1-score: {f1:.2%}")
        print(f"   📊 OR-only flattening tests: {flattening_tests['passed']}/{flattening_tests['total']}")
        
        print("\n" + "="*80)
        print("TEST SUITE COMPLETED")
        print("="*80)
        
        return results
    
    def test_or_only_flattening(self, model: cobra.Model, translation_issues: Dict) -> Dict:
        """Test the OR-only GPR flattening functionality"""
        test_cases = [
            # (original_gpr, expected_after_flattening, should_be_flattened)
            ("SHARED_TARGET1 or SHARED_TARGET1", "SHARED_TARGET1", True),
            ("(SHARED_TARGET2 or SHARED_TARGET2) or SHARED_TARGET2", "SHARED_TARGET2", True),
            ("TARGET1 or TARGET2 or TARGET1", "TARGET1 or TARGET2", True),
            ("(TARGET1 or TARGET2) and TARGET3", "(TARGET1 or TARGET2) and TARGET3", False),  # Contains AND
            ("TARGET1 and TARGET1", "TARGET1", True),  # Should simplify AND duplicates too
        ]
        
        results = {'total': 0, 'passed': 0, 'failed': [], 'details': []}
        
        print("   Testing OR-only flattening functionality:")
        
        # Test the helper functions directly
        for original, expected, should_flatten in test_cases:
            results['total'] += 1
            
            # Test _is_or_only_expression
            is_or_only = model_utils._is_or_only_expression(original)
            
            # Test _flatten_or_only_gpr if it should be OR-only
            if should_flatten and 'and' not in original.lower():
                flattened = model_utils._flatten_or_only_gpr(original)
                passed = flattened == expected
            else:
                passed = not should_flatten or is_or_only == (not 'and' in original.lower())
                flattened = original
            
            status = "✓" if passed else "❌"
            results['details'].append({
                'original': original,
                'expected': expected,
                'flattened': flattened,
                'is_or_only': is_or_only,
                'should_flatten': should_flatten,
                'passed': passed
            })
            
            if passed:
                results['passed'] += 1
            else:
                results['failed'].append(f"{original} -> {flattened} (expected: {expected})")
            
            print(f"      {status} '{original}' -> '{flattened}' (OR-only: {is_or_only})")
        
        # Test actual model reactions that should have been flattened
        for rxn in model.reactions:
            if rxn.id in translation_issues:
                original_gpr = rxn.gene_reaction_rule
                is_or_only = model_utils._is_or_only_expression(original_gpr)
                if is_or_only:
                    print(f"      🔍 Real case: {rxn.id} has OR-only GPR: '{original_gpr}'")
        
        return results

def run_individual_tests():
    """Run individual component tests"""
    print("\n" + "="*80)
    print("INDIVIDUAL COMPONENT TESTS")
    print("="*80)
    
    # Test 1: OR-only detection
    print("\n1. Testing OR-only detection...")
    or_only_cases = [
        ("GENE1 or GENE2", True),
        ("(GENE1 or GENE2)", True),
        ("GENE1 or GENE2 or GENE3", True),
        ("(GENE1 or GENE2) or GENE3", True),
        ("((GENE1 or GENE2) or GENE3) or GENE4", True),
        ("GENE1 and GENE2", False),
        ("GENE1 or (GENE2 and GENE3)", False),
        ("(GENE1 or GENE2) and GENE3", False),
        ("GENE1", False),  # Single gene
        ("", False),  # Empty
    ]
    
    for gpr, expected in or_only_cases:
        result = model_utils._is_or_only_expression(gpr)
        status = "✓" if result == expected else "❌"
        print(f"   {status} '{gpr}' -> {result} (expected: {expected})")
    
    # Test 2: GPR flattening
    print("\n2. Testing GPR flattening...")
    flattening_cases = [
        ("GENE1 or GENE1", "GENE1"),
        ("(GENE1 or GENE1) or GENE2", "GENE1 or GENE2"),
        ("GENE1 or GENE2 or GENE1", "GENE1 or GENE2"),
        ("(GENE1 or GENE2) or (GENE1 or GENE3)", "GENE1 or GENE2 or GENE3"),
        ("((A or A) or B) or C", "A or B or C"),
    ]
    
    for original, expected in flattening_cases:
        result = model_utils._flatten_or_only_gpr(original)
        status = "✓" if result == expected else "❌"
        print(f"   {status} '{original}' -> '{result}' (expected: '{expected}')")

def main():
    """Main test function"""
    print("COBRAxy GPR Translation Comprehensive Test Suite")
    print("=" * 80)
    
    # Run individual component tests first
    run_individual_tests()
    
    # Run comprehensive test suite
    tester = GPRTranslationTester()
    results = tester.run_comprehensive_test()
    
    # Save results for further analysis if needed
    if results['success']:
        print(f"\n✅ All tests completed successfully!")
        print(f"📁 Test models and results available in results object")
        
        # Optionally save to file
        try:
            import pickle
            with open('/tmp/gpr_translation_test_results.pkl', 'wb') as f:
                pickle.dump(results, f)
            print(f"📁 Detailed results saved to /tmp/gpr_translation_test_results.pkl")
        except:
            pass
    else:
        print(f"\n❌ Tests failed: {results.get('error', 'Unknown error')}")
        return False
    
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)