#pragma once

namespace plksim {
namespace io {

enum OutputFormat {
  dx,      // OpenDX
  eps,     // Encapsulated PostScript
  gmv,     // General Mesh Viewer
  gnuplot, // Gnuplot
  hdf5,    // Hierarchical Data Format
  povray,  // Persistence of Vision Raytracer
  svg,     // SVG
  tecplot, // Tecplot
  ucd,     // UCD for AVS
  vtk,     // VTK
  vtu      // VTU
};

}
} // namespace plksim