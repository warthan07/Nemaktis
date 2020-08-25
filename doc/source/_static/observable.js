import {Runtime,Library} from "https://cdn.jsdelivr.net/npm/@observablehq/runtime@4/dist/runtime.js";

export function getRuntime(fig_id) {
    const stdlib = new Library;
    const target = document.querySelector(fig_id);

    function width() {
		return stdlib.Generators.observe(notify => {
			let width = notify(target.clientWidth);
			function resized() {
				let width1 = target.clientWidth;
				if (width1 !== width) notify(width = width1);
			}
			window.addEventListener("resize", resized);
			return () => window.removeEventListener("resize", resized);
		});
    }

	return (new Runtime(Object.assign(stdlib, {width:width})));
}
