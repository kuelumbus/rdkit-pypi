"""
Test suite for RDKit Python stub files (.pyi).
"""

from pathlib import Path
import pytest


class TestRDKitStubs:
    """Test RDKit Python stub files and type hints."""

    def test_stub_files_exist(self):
        """Test that key .pyi stub files exist in the installation."""
        # Try to find stub files in the installed package
        try:
            import rdkit
            rdkit_path = Path(rdkit.__file__).parent.parent
            stubs_path = rdkit_path / "rdkit-stubs"
            
            if not stubs_path.exists():
                pytest.fail("rdkit-stubs directory not found in installation")
            
            # Check for key stub files
            expected_stubs = [
                "Chem/__init__.pyi",
                "Chem/Descriptors.pyi", 
                "Chem/rdMolDescriptors.pyi",
                "Chem/AllChem.pyi",
            ]
            
            missing_stubs = []
            for stub_file in expected_stubs:
                stub_path = stubs_path / stub_file
                if not stub_path.exists():
                    missing_stubs.append(stub_file)
            
            if missing_stubs:
                pytest.fail(f"Missing stub files: {missing_stubs}")
                
        except ImportError:
            pytest.fail("rdkit package not available for stub testing")






