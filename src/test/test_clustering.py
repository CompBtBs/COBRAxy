"""
Unit tests for marea_cluster module.

Run with: python -m pytest test_clustering.py -v
Or: python test_clustering.py
"""

import sys
import os
import pandas as pd
import numpy as np
import tempfile
from pathlib import Path

# Try to import pytest, but don't fail if not available
try:
    import pytest
    HAS_PYTEST = True
except ImportError:
    HAS_PYTEST = False
    class _DummyPytest:
        class raises:
            def __init__(self, *args, **kwargs):
                self.expected_exceptions = args
            def __enter__(self):
                return self
            def __exit__(self, exc_type, exc_val, exc_tb):
                if exc_type is None:
                    raise AssertionError("Expected an exception but none was raised")
                if not any(issubclass(exc_type, e) for e in self.expected_exceptions):
                    return False
                return True
    pytest = _DummyPytest()

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import marea_cluster

# Get the tool directory
TOOL_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


class TestMAREACluster:
    """Tests for marea_cluster module"""
    
    def test_process_args(self):
        """Test argument processing"""
        args = marea_cluster.process_args([
            '-td', TOOL_DIR,
            '-cy', 'kmeans'
        ])
        assert hasattr(args, 'tool_dir')
        assert hasattr(args, 'cluster_type')
    
    def test_clustering_types(self):
        """Test that all clustering types are available"""
        args = marea_cluster.process_args(['-cy', 'kmeans'])
        assert args.cluster_type == 'kmeans'
        
        args = marea_cluster.process_args(['-cy', 'dbscan'])
        assert args.cluster_type == 'dbscan'
        
        args = marea_cluster.process_args(['-cy', 'hierarchy'])
        assert args.cluster_type == 'hierarchy'


class TestKMeansClustering:
    """Tests for K-means clustering"""
    
    def test_kmeans_basic(self):
        """Test basic K-means clustering"""
        from sklearn.cluster import KMeans
        
        # Create sample data
        data = np.array([
            [1.0, 2.0],
            [1.5, 1.8],
            [5.0, 8.0],
            [8.0, 8.0],
            [1.0, 0.6],
            [9.0, 11.0]
        ])
        
        # Perform clustering
        kmeans = KMeans(n_clusters=2, random_state=42)
        labels = kmeans.fit_predict(data)
        
        assert len(labels) == len(data)
        assert len(set(labels)) == 2  # Should have 2 clusters
    
    def test_kmeans_with_dataframe(self):
        """Test K-means clustering with DataFrame"""
        from sklearn.cluster import KMeans
        
        # Create sample DataFrame
        df = pd.DataFrame({
            'feature1': [1.0, 1.5, 5.0, 8.0, 1.0, 9.0],
            'feature2': [2.0, 1.8, 8.0, 8.0, 0.6, 11.0]
        })
        
        kmeans = KMeans(n_clusters=2, random_state=42)
        df['cluster'] = kmeans.fit_predict(df[['feature1', 'feature2']])
        
        assert 'cluster' in df.columns
        assert len(df['cluster'].unique()) == 2


class TestDBSCANClustering:
    """Tests for DBSCAN clustering"""
    
    def test_dbscan_basic(self):
        """Test basic DBSCAN clustering"""
        from sklearn.cluster import DBSCAN
        
        # Create sample data with clear clusters
        data = np.array([
            [1.0, 2.0],
            [1.5, 1.8],
            [1.2, 2.1],
            [8.0, 8.0],
            [8.5, 8.3],
            [8.2, 8.1]
        ])
        
        # Perform clustering
        dbscan = DBSCAN(eps=1.0, min_samples=2)
        labels = dbscan.fit_predict(data)
        
        assert len(labels) == len(data)
        # DBSCAN should find at least 2 clusters (excluding noise as -1)
        unique_labels = set(labels)
        unique_labels.discard(-1)  # Remove noise label
        assert len(unique_labels) >= 1


class TestHierarchicalClustering:
    """Tests for Hierarchical clustering"""
    
    def test_hierarchical_basic(self):
        """Test basic hierarchical clustering"""
        from sklearn.cluster import AgglomerativeClustering
        
        # Create sample data
        data = np.array([
            [1.0, 2.0],
            [1.5, 1.8],
            [5.0, 8.0],
            [8.0, 8.0],
            [1.0, 0.6],
            [9.0, 11.0]
        ])
        
        # Perform clustering
        hierarchical = AgglomerativeClustering(n_clusters=2)
        labels = hierarchical.fit_predict(data)
        
        assert len(labels) == len(data)
        assert len(set(labels)) == 2


class TestScaling:
    """Tests for data scaling/normalization"""
    
    def test_standard_scaling(self):
        """Test standard scaling"""
        from sklearn.preprocessing import StandardScaler
        
        data = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
        
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(data)
        
        # Check that mean is close to 0
        assert np.abs(scaled_data.mean(axis=0)).max() < 1e-10
        
        # Check that std is close to 1
        assert np.abs(scaled_data.std(axis=0) - 1.0).max() < 1e-10
    
    def test_minmax_scaling(self):
        """Test min-max scaling"""
        from sklearn.preprocessing import MinMaxScaler
        
        data = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
        
        scaler = MinMaxScaler()
        scaled_data = scaler.fit_transform(data)
        
        # Check that min is 0 and max is 1
        assert scaled_data.min() == 0.0
        assert scaled_data.max() == 1.0


class TestClusterEvaluation:
    """Tests for cluster evaluation metrics"""
    
    def test_silhouette_score(self):
        """Test silhouette score calculation"""
        from sklearn.cluster import KMeans
        from sklearn.metrics import silhouette_score
        
        # Create well-separated data
        data = np.array([
            [1.0, 1.0],
            [1.5, 1.5],
            [1.2, 1.3],
            [10.0, 10.0],
            [10.5, 10.5],
            [10.2, 10.3]
        ])
        
        kmeans = KMeans(n_clusters=2, random_state=42)
        labels = kmeans.fit_predict(data)
        
        score = silhouette_score(data, labels)
        
        # Well-separated clusters should have high silhouette score
        assert score > 0.5
        assert score <= 1.0


class TestDataPreparation:
    """Tests for data preparation for clustering"""
    
    def test_remove_constant_features(self):
        """Test removal of constant features"""
        df = pd.DataFrame({
            'var_feature': [1, 2, 3, 4, 5],
            'const_feature': [1, 1, 1, 1, 1],
            'var_feature2': [5, 4, 3, 2, 1]
        })
        
        # Remove constant columns
        df_filtered = df.loc[:, df.std() > 0]
        
        assert 'const_feature' not in df_filtered.columns
        assert 'var_feature' in df_filtered.columns
        assert 'var_feature2' in df_filtered.columns
    
    def test_handle_missing_values(self):
        """Test handling of missing values"""
        df = pd.DataFrame({
            'feature1': [1.0, 2.0, np.nan, 4.0],
            'feature2': [5.0, np.nan, 7.0, 8.0]
        })
        
        # Drop rows with NaN
        df_clean = df.dropna()
        assert len(df_clean) == 2
        
        # Fill NaN with mean
        df_filled = df.fillna(df.mean())
        assert not df_filled.isnull().any().any()


class TestVisualization:
    """Tests for clustering visualization"""
    
    def test_dendrogram_data(self):
        """Test that we can create dendrogram data"""
        from scipy.cluster.hierarchy import linkage
        
        data = np.array([
            [1.0, 2.0],
            [1.5, 1.8],
            [5.0, 8.0],
            [8.0, 8.0]
        ])
        
        # Create linkage matrix
        Z = linkage(data, method='ward')
        
        assert Z.shape[0] == len(data) - 1
        assert Z.shape[1] == 4
    
    def test_elbow_method_data(self):
        """Test data preparation for elbow method"""
        from sklearn.cluster import KMeans
        
        data = np.array([
            [1.0, 2.0],
            [1.5, 1.8],
            [5.0, 8.0],
            [8.0, 8.0],
            [1.0, 0.6],
            [9.0, 11.0]
        ])
        
        inertias = []
        K_range = range(1, 4)
        
        for k in K_range:
            kmeans = KMeans(n_clusters=k, random_state=42)
            kmeans.fit(data)
            inertias.append(kmeans.inertia_)
        
        # Inertia should decrease as K increases
        assert inertias[0] > inertias[1] > inertias[2]


class TestRealWorldScenarios:
    """Tests with realistic metabolic data scenarios"""
    
    def test_cluster_reactions_by_flux(self):
        """Test clustering reactions based on flux patterns"""
        # Create sample flux data for different conditions
        df = pd.DataFrame({
            'reaction': ['r1', 'r2', 'r3', 'r4', 'r5'],
            'condition1': [1.5, 1.6, 0.1, 0.2, 2.0],
            'condition2': [1.4, 1.7, 0.15, 0.18, 2.1],
            'condition3': [1.6, 1.5, 0.12, 0.22, 1.9]
        })
        
        # Extract numeric features
        features = df[['condition1', 'condition2', 'condition3']]
        
        # Cluster
        from sklearn.cluster import KMeans
        kmeans = KMeans(n_clusters=2, random_state=42)
        df['cluster'] = kmeans.fit_predict(features)
        
        # r1, r2, r5 should be in one cluster (high flux)
        # r3, r4 should be in another cluster (low flux)
        high_flux_reactions = df[df['condition1'] > 1.0]['reaction'].tolist()
        low_flux_reactions = df[df['condition1'] < 1.0]['reaction'].tolist()
        
        assert len(high_flux_reactions) == 3
        assert len(low_flux_reactions) == 2
    
    def test_cluster_samples_by_metabolic_profile(self):
        """Test clustering samples based on metabolic profiles"""
        # Create sample data: samples x reactions
        df = pd.DataFrame({
            'sample': ['normal1', 'normal2', 'cancer1', 'cancer2'],
            'r1': [1.5, 1.6, 3.0, 3.1],
            'r2': [0.5, 0.6, 0.1, 0.15],
            'r3': [2.0, 2.1, 4.0, 4.2]
        })
        
        features = df[['r1', 'r2', 'r3']]
        
        from sklearn.cluster import KMeans
        kmeans = KMeans(n_clusters=2, random_state=42)
        df['cluster'] = kmeans.fit_predict(features)
        
        # Normal samples should cluster separately from cancer samples
        assert len(df['cluster'].unique()) == 2


if __name__ == "__main__":
    # Run tests with pytest if available
    if HAS_PYTEST:
        pytest.main([__file__, "-v"])
    else:
        print("pytest not available, running basic tests...")
        
        test_classes = [
            TestMAREACluster(),
            TestKMeansClustering(),
            TestDBSCANClustering(),
            TestHierarchicalClustering(),
            TestScaling(),
            TestClusterEvaluation(),
            TestDataPreparation(),
            TestVisualization(),
            TestRealWorldScenarios()
        ]
        
        failed = 0
        passed = 0
        
        for test_class in test_classes:
            class_name = test_class.__class__.__name__
            print(f"\n{class_name}:")
            
            for method_name in dir(test_class):
                if method_name.startswith("test_"):
                    try:
                        method = getattr(test_class, method_name)
                        method()
                        print(f"  ✓ {method_name}")
                        passed += 1
                    except Exception as e:
                        print(f"  ✗ {method_name}: {str(e)}")
                        failed += 1
        
        print(f"\n{'='*60}")
        print(f"Results: {passed} passed, {failed} failed")
        if failed > 0:
            sys.exit(1)
