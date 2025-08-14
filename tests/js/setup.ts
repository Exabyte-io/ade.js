import JSONSchemasInterface from "@mat3ra/esse/dist/js/esse/JSONSchemasInterface";
import type { JSONSchema } from "@mat3ra/esse/dist/js/esse/utils";
import schemas from "@mat3ra/esse/dist/js/schemas.json";

// Global setup that runs once before all tests
before(() => {
    JSONSchemasInterface.setSchemas(schemas as JSONSchema[]);
});
