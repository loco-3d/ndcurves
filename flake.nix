{
  description = "Library for creating smooth cubic splines";

  inputs = {
    gepetto.url = "github:gepetto/nix";
    flake-parts.follows = "gepetto/flake-parts";
    systems.follows = "gepetto/systems";
  };

  outputs =
    inputs:
    inputs.flake-parts.lib.mkFlake { inherit inputs; } (
      { lib, ... }:
      {
        systems = import inputs.systems;
        imports = [
          inputs.gepetto.flakeModule
          {
            flakoboros = {
              extraDevPyPackages = [ "ndcurves" ];
              overrideAttrs.ndcurves = _: {
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
            };
          }
        ];
      }
    );
}
