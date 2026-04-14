{
  description = "Library for creating smooth cubic splines";

  inputs.gepetto.url = "github:gepetto/nix";

  outputs =
    inputs:
    inputs.gepetto.lib.mkFlakoboros inputs (
      { lib, ... }:
      {
        extraDevPyPackages = [ "ndcurves" ];
        overrideAttrs.ndcurves = {
          src = lib.fileset.toSource {
            root = ./.;
            fileset = lib.fileset.unions [
              ./CMakeLists.txt
              ./doc
              ./include
              ./package.xml
              ./python
              ./tests
            ];
          };
        };
      }
    );
}
