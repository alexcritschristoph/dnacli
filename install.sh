#!/usr/bin/env sh
# Installs dnacli from its GitHub releases.
#
#   curl -fsSL https://raw.githubusercontent.com/alexcritschristoph/dnacli/main/install.sh | sh
#
# DNACLI_BIN_DIR=/usr/local/bin puts dnacli somewhere other than ~/.local/bin.
# DNACLI_VERSION=v1.2.3 installs that release instead of the latest one.
set -eu

repo=alexcritschristoph/dnacli
version=${DNACLI_VERSION:-latest}
bindir=${DNACLI_BIN_DIR:-$HOME/.local/bin}

case $(uname -s) in
Linux) goos=linux ;;
Darwin) goos=darwin ;;
*)
	echo "dnacli runs on Linux and macOS. $(uname -s) is neither." >&2
	exit 1
	;;
esac

case $(uname -m) in
x86_64 | amd64) goarch=amd64 ;;
arm64 | aarch64) goarch=arm64 ;;
*)
	echo "dnacli is built for x86-64 and ARM64. $(uname -m) is neither." >&2
	exit 1
	;;
esac

if [ "$version" = latest ]; then
	base=https://github.com/$repo/releases/latest/download
else
	base=https://github.com/$repo/releases/download/$version
fi
name=dnacli-$goos-$goarch

if command -v curl >/dev/null; then
	download() { curl -fsSL "$1"; }
else
	download() { wget -qO- "$1"; }
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

echo "Downloading $name.tar.gz."
if ! download "$base/$name.tar.gz" > "$tmp/$name.tar.gz"; then
	echo "Could not download $base/$name.tar.gz" >&2
	exit 1
fi

want=$(download "$base/checksums.txt" | awk -v f="$name.tar.gz" '$2 == f { print $1 }')
if [ -z "$want" ]; then
	echo "The release has no checksum for $name.tar.gz. Not installing." >&2
	exit 1
fi
if command -v sha256sum >/dev/null; then
	got=$(sha256sum "$tmp/$name.tar.gz" | awk '{ print $1 }')
else
	got=$(shasum -a 256 "$tmp/$name.tar.gz" | awk '{ print $1 }')
fi
if [ "$got" != "$want" ]; then
	echo "The download is damaged; its checksum does not match. Not installing." >&2
	exit 1
fi

tar -xzf "$tmp/$name.tar.gz" -C "$tmp"
mkdir -p "$bindir"
mv "$tmp/$name/dnacli" "$bindir/dnacli"
chmod +x "$bindir/dnacli"
echo "Installed $bindir/dnacli"

case ":$PATH:" in
*":$bindir:"*) ;;
*)
	echo
	echo "$bindir is not in your PATH. Add it to your shell's startup file:"
	echo
	echo "    export PATH=\"\$PATH:$bindir\""
	;;
esac
