{
  description = "PotreeConverter - multi-res point cloud converter";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};

        PotreeConverter = pkgs.callPackage ./package.nix { };
      in
      {
        # nix build
        packages = {
          inherit PotreeConverter;
          default = PotreeConverter;
        };

        # nix develop
        devShells.default = pkgs.mkShell {
          name = "PotreeConverter-dev";

          inputsFrom = [ PotreeConverter ];

          packages = with pkgs; [
            cmake
            ninja
            clang-tools 
            gdb
          ];

          CMAKE_POLICY_VERSION_MINIMUM="3.5";
        };
      }
    );
}
