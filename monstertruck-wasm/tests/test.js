import init, * as Truck from "../pkg/truck_js.js";
import { assertEquals } from "https://deno.land/std@0.110.0/testing/asserts.ts";

function assertNearArray(actual, expected, tolerance = 1.0e-6) {
  assertEquals(actual.length, expected.length);
  for (let i = 0; i < actual.length; i++) {
    if (Math.abs(actual[i] - expected[i]) > tolerance) {
      throw new Error(
        `component ${i}: got ${actual[i]}, expected ${expected[i]} within ${tolerance}`,
      );
    }
  }
}

function roundTripSolid(solid) {
  const decoded = Truck.Solid.from_json(solid.to_json());
  if (decoded === undefined || decoded === null) {
    throw new Error("Solid.from_json must accept Solid.to_json output");
  }
  return decoded;
}

Deno.test("swept cube has unit-cube geometry and JSON round-trips", async () => {
  await init(Deno.readFile("./pkg/truck_js_bg.wasm"));
  const v = Truck.vertex(-0.5, -0.5, -0.5);
  const e = Truck.extrude(v.upcast(), [1.0, 0.0, 0.0]);
  const f = Truck.extrude(e, [0.0, 1.0, 0.0]);
  const abst = Truck.extrude(f, [0.0, 0.0, 1.0]);
  const solid = abst.into_solid();

  const bbox = solid.to_polygon(0.01).bounding_box();
  assertNearArray(bbox, [-0.5, -0.5, -0.5, 0.5, 0.5, 0.5]);

  const decoded = roundTripSolid(solid);
  assertNearArray(decoded.to_polygon(0.01).bounding_box(), bbox);
});

Deno.test("revolved torus JSON round-trip preserves geometry", async () => {
  await init(Deno.readFile("./pkg/truck_js_bg.wasm"));
  const v = Truck.vertex(0.75, 0.0, 0.25);
  const w = Truck.revolve(v.upcast(), [0.75, 0.0, 0.0], [0.0, 1.0, 0.0], 7.0, 2);
  const abst = Truck.revolve(w, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 7.0, 2);
  const shell = abst.into_shell();
  const solid = shell.into_solid();

  const bbox = solid.to_polygon(0.01).bounding_box();
  const decoded = roundTripSolid(solid);
  assertNearArray(decoded.to_polygon(0.01).bounding_box(), bbox);
});
