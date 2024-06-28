{
  lib,
  cmake,
  jrl-cmakemodules,
  python3Packages,
}:

python3Packages.buildPythonPackage {
  pname = "ndcurves";
  version = "1.4.1";
  pyproject = false;

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

  strictDeps = true;

  nativeBuildInputs = [ cmake ];
  propagatedBuildInputs = [
    jrl-cmakemodules
    python3Packages.eigenpy
    python3Packages.pinocchio
  ];

  cmakeFlags = [ "-DCURVES_WITH_PINOCCHIO_SUPPORT=ON" ];

  doCheck = true;

  pythonImportsCheck = [ "ndcurves" ];

  meta = {
    description = "Environments and robot descriptions for HPP";
    homepage = "https://github.com/humanoid-path-planner/hpp-environments";
    license = lib.licenses.bsd2;
    maintainers = [ lib.maintainers.nim65s ];
  };
}
