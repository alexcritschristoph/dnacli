#!/usr/bin/env bash
# Builds dnacli for every computer it ships to.
#
# Each target gets its own folder holding a binary named dnacli, so nobody has
# to rename what they download, and a tarball of that folder goes in dist/.
# ./build.sh 1.2.3 stamps that version into the binary. Without an argument the
# version comes from a git tag, and failing that the one in main.go stands.
set -euo pipefail

version=${1:-$(git describe --tags --dirty 2>/dev/null || true)}
version=${version#v}
dist=dist

ldflags="-s -w"
if [ -n "$version" ]; then
	ldflags="$ldflags -X main.version=$version"
fi

rm -rf "$dist"
mkdir "$dist"

for target in linux/amd64 linux/arm64 darwin/amd64 darwin/arm64; do
	goos=${target%/*}
	goarch=${target#*/}
	name=dnacli-$goos-$goarch
	CGO_ENABLED=0 GOOS=$goos GOARCH=$goarch \
		go build -trimpath -ldflags "$ldflags" -o "$dist/$name/dnacli" ./src
	tar -czf "$dist/$name.tar.gz" -C "$dist" "$name"
	rm -r "$dist/$name"
done

cd "$dist"
if command -v sha256sum >/dev/null; then
	sha256sum *.tar.gz > checksums.txt
else
	shasum -a 256 *.tar.gz > checksums.txt
fi

echo "Built:"
ls -1
