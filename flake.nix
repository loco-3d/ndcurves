{
  description = "Library for creating smooth cubic splines";

  inputs = {
    gepetto.url = "github:gepetto/nix";
    flake-parts.follows = "gepetto/flake-parts";
    nixpkgs.follows = "gepetto/nixpkgs";
    nix-ros-overlay.follows = "gepetto/nix-ros-overlay";
    systems.follows = "gepetto/systems";
    treefmt-nix.follows = "gepetto/treefmt-nix";
  };

  outputs =
    inputs:
    inputs.flake-parts.lib.mkFlake { inherit inputs; } {
      systems = import inputs.systems;
      imports = [ inputs.gepetto.flakeModule ];
      perSystem =
        {
          lib,
          pkgs,
          self',
          ...
        }:
        {
          apps.default = {
            type = "app";
            program = pkgs.python3.withPackages (_: [ self'.packages.default ]);
          };
          packages = {
            default = self'.packages.py-ndcurves;
            ndcurves = pkgs.ndcurves.overrideAttrs {
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
            py-ndcurves = pkgs.python3Packages.toPythonModule (
              self'.packages.ndcurves.overrideAttrs (super: {
                pname = "py-${super.pname}";
                cmakeFlags = super.cmakeFlags ++ [
                  (lib.cmakeBool "BUILD_PYTHON_INTERFACE" true)
                  (lib.cmakeBool "BUILD_STANDALONE_PYTHON_INTERFACE" true)
                ];
                nativeBuildInputs = super.nativeBuildInputs ++ [
                  pkgs.python3Packages.python
                ];
                propagatedBuildInputs = [
                  pkgs.python3Packages.pinocchio
                  self'.packages.ndcurves
                ] ++ super.propagatedBuildInputs;
                nativeCheckInputs = [
                  pkgs.python3Packages.pythonImportsCheckHook
                ];
              })
            );
          };
        };
    };
}
