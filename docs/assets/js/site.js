// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

// Importing JavaScript
// Import Bootstrap's bundle (all of Bootstrap's JS + Popper.js dependency)
import "../bootstrap/bootstrap.bundle.min.js";

async function loadIconSprite(path = "./assets/icons.svg") {
  const res = await fetch(path);
  if (!res.ok) throw new Error(`Icon sprite load failed: ${res.status} ${res.statusText}`);

  const svgText = await res.text();
  const wrapper = document.createElement("div");
  wrapper.innerHTML = svgText;

  const svg = wrapper.querySelector("svg");
  if (!svg) throw new Error("icons.svg did not contain an <svg> root");

  svg.style.display = "none";
  document.body.prepend(svg);
}

loadIconSprite().catch(console.error);
