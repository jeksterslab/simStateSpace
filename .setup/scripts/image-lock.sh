#!/bin/bash
# The image job of the GitHub Actions workflows (rProject 0.0.57): reads
# .setup/scripts/rproject.lock (written by make auto) and outputs the image
# to run (`image`), the Posit Package Manager snapshot (`snapshot`), and the
# R version (`r_version`). A scheduled run, a manual run with `latest`, no
# lock, or a lock without an image run the latest image and packages.
#   EVENT   github.event_name
#   LATEST  the workflow_dispatch input `latest` (true or false)
#   IMAGE   the latest image (default ijapesigan/docs)

set -u

lock=".setup/scripts/rproject.lock"
latest="${IMAGE:-ijapesigan/docs}"
image="${latest}"
snapshot="latest"
r_version="release"

field() {
  sed -n "s/^$1 <- \"\(.*\)\"$/\1/p" "${lock}" 2> /dev/null | head -n 1
}

if [ "${EVENT:-}" = "schedule" ]; then
  echo "::notice::A scheduled run uses the latest image (${latest})."
elif [ "${LATEST:-false}" = "true" ]; then
  echo "::notice::latest requested: running ${latest}."
elif [ ! -f "${lock}" ]; then
  echo "::notice::No ${lock}: running ${latest}."
else
  locked="$(field image)"
  name="$(field image_name)"
  if [ -z "${locked}" ]; then
    echo "::notice::The lock records no image: running ${latest}."
  elif [[ "${locked}" == *@sha256:* ]]; then
    if docker manifest inspect "${locked}" > /dev/null 2>&1; then
      image="${locked}"
      snapshot="$(field snapshot)"
      r_version="$(field r_version)"
    else
      image="${name:-${latest}}"
      echo "::warning::${locked} is gone from the registry: running ${image}:latest."
    fi
  else
    image="${locked}"
    snapshot="$(field snapshot)"
    r_version="$(field r_version)"
    echo "::warning::The lock has no image digest: running ${image}:latest. Rebuild the add-on sif from the public sif (Bootstrap: localimage)."
  fi
fi

{
  echo "image=${image}"
  echo "snapshot=${snapshot:-latest}"
  echo "r_version=${r_version:-release}"
} >> "${GITHUB_OUTPUT}"
echo "Image: ${image}; packages: ${snapshot:-latest}; R: ${r_version:-release}"
