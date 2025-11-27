#!/usr/bin/env python3
"""
Element Implementation Verification Script
============================================
Verifies that all element classes have both header and implementation files.
"""

import os
from pathlib import Path

# Element definitions: (class_name, element_type)
ELEMENTS = [
    # 3D Solid Elements (Linear)
    ("RgTet4Element", "Tet4 - 4-node linear tetrahedral"),
    ("RgHex8Element", "Hex8 - 8-node linear hexahedral"),
    
    # 3D Solid Elements (Quadratic)
    ("RgTet10Element", "Tet10 - 10-node quadratic tetrahedral"),
    ("RgHex20Element", "Hex20 - 20-node serendipity hexahedral"),
    
    # 3D Nonlinear Elements
    ("RgHex8GeomNLElement", "Hex8NL - 8-node geometric nonlinear hex"),
    
    # 1D Structural Elements
    ("RgBeam2dElement", "Beam2D - 2-node 2D Timoshenko beam"),
    ("RgBeam3dElement", "Beam3D - 2-node 3D Timoshenko beam"),
]

def check_element_files(base_path="e:\\Lzs\\FEM\\lzsfem\\RgFem\\src\\elements\\RgElement"):
    """Check if all element files exist."""
    
    print("=" * 80)
    print("ELEMENT IMPLEMENTATION VERIFICATION")
    print("=" * 80)
    print()
    
    all_complete = True
    completed = []
    missing = []
    
    for class_name, description in ELEMENTS:
        header_file = Path(base_path) / f"{class_name}.h"
        impl_file = Path(base_path) / f"{class_name}.cpp"
        
        header_exists = header_file.exists()
        impl_exists = impl_file.exists()
        
        status = "✅" if (header_exists and impl_exists) else "❌"
        
        print(f"{status} {class_name}")
        print(f"   Type: {description}")
        print(f"   Header: {'✓' if header_exists else '✗'} {class_name}.h")
        print(f"   Implementation: {'✓' if impl_exists else '✗'} {class_name}.cpp")
        
        if header_exists and impl_exists:
            completed.append(class_name)
            # Get file sizes
            header_size = header_file.stat().st_size
            impl_size = impl_file.stat().st_size
            print(f"   Size: H({header_size:,}B) + C({impl_size:,}B)")
        else:
            all_complete = False
            missing.append((class_name, not header_exists, not impl_exists))
        
        print()
    
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total Elements: {len(ELEMENTS)}")
    print(f"Completed: {len(completed)}/{len(ELEMENTS)}")
    
    if completed:
        print("\n✅ COMPLETED IMPLEMENTATIONS:")
        for elem in completed:
            print(f"   - {elem}")
    
    if missing:
        print("\n❌ MISSING IMPLEMENTATIONS:")
        for elem, no_header, no_impl in missing:
            missing_parts = []
            if no_header:
                missing_parts.append("header")
            if no_impl:
                missing_parts.append("implementation")
            print(f"   - {elem} (missing {', '.join(missing_parts)})")
    
    print("\n" + "=" * 80)
    if all_complete:
        print("✅ ALL ELEMENT CLASSES ARE FULLY IMPLEMENTED!")
        print("=" * 80)
        return True
    else:
        print("⚠️  Some implementations are still pending.")
        print("=" * 80)
        return False

if __name__ == "__main__":
    success = check_element_files()
    exit(0 if success else 1)
