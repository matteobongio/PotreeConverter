# https://github.com/NixOS/nixpkgs/blob/nixos-25.11/pkgs/by-name/po/potreeconverter/package.nix
{
  lib,
  stdenv,
  cmake,
  boost,
  onetbb,
  makeWrapper,
}:

stdenv.mkDerivation {
  pname = "PotreeConverter";
  version = "unstable-2023-02-27";

  src = ./.;

  buildInputs = [
    boost
    onetbb
  ];

  nativeBuildInputs = [
    makeWrapper
    cmake
  ];

  cmakeFlags = [ (lib.strings.cmakeFeature "CMAKE_POLICY_VERSION_MINIMUM" "3.5") ];

  postPatch = ''
    runHook prePatch

    substituteInPlace ./CMakeLists.txt \
      --replace "find_package(TBB REQUIRED)" ""

    # prevent inheriting permissions from /nix/store when copying
    substituteInPlace Converter/src/main.cpp --replace \
      'fs::copy(templateDir, pagedir, fs::copy_options::overwrite_existing | fs::copy_options::recursive)' 'string cmd = "cp --no-preserve=mode -r " + templateDir + " " + pagedir; system(cmd.c_str());'
  '';

  # The upstream build system does not provide an install target.
  installPhase = ''
    runHook preInstall

    mkdir -p $out/{bin,lib}
    mv liblaszip.so $out/lib
    mv PotreeConverter $out/bin
    ln -s $out/bin/PotreeConverter $out/bin/potreeconverter

    # Create an empty wrapper, since PotreeConverter segfaults if called via
    # $PATH rather than absolute path. An empty wrapper forces an absolute path
    # on each invocation
    wrapProgram $out/bin/PotreeConverter

    runHook postInstall
  '';

  postFixup = ''
    ln -s $src/resources $out/bin/resources
  '';

  meta = {
    description = "Create multi res point cloud to use with potree";
    homepage = "https://github.com/potree/PotreeConverter";
    license = lib.licenses.bsd2;
    maintainers = with lib.maintainers; [ matthewcroughan ];
    platforms = with lib.platforms; linux;
  };
}
