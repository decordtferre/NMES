// You should not need to modify this file

// This Jacobian refers to the transformation between
// the reference element (hard coded in the finite-element code)
// and the element actual in the mesh

// Note that the default Jacobian is Vol and corresponds to the higher dimension.
// In a 3D model, there is no confusion: the higher dimension is N=3, geometry is split in
// volumes (tetrahedra, hexahedra, ...), and the Jacobian is Vol.
// If we have one equation explicitly defined on a boundary, i.e. a surface, the
// Jacobian there is associated to the dimension N-1=2, surfaces, it is Sur
//
// In a 2D model, the higher dimension is N=2 (discretisation in triangles, quadrangles...),
// however the Jacobian associated to surfaces is Vol,
// the one associated to the dimension N-1=1, lines is Sur.

Jacobian {
  { Name Vol;
    Case {
      { Region All; Jacobian Vol; }
    }
  }
  { Name Sur;
    Case {
      { Region All; Jacobian Sur; }
    }
  }
}


// What follows refer to the number of Gauss integration points to evaluate the integrals per finite element
// Flexibility is allowed for research reasons. The defaults hereafter should be sufficient for you
Integration {
  { Name I1;
    Case {
      { Type Gauss;
        Case {
          { GeoElement Point; NumberOfPoints  1; }
          { GeoElement Line; NumberOfPoints  4; }
          { GeoElement Triangle; NumberOfPoints  4; }
          { GeoElement Quadrangle; NumberOfPoints  4; }

          { GeoElement Triangle2; NumberOfPoints  7; } // Second order, geometrical elements
        }
      }
    }
  }
}
