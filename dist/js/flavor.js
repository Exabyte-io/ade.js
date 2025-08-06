"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const entity_1 = require("@mat3ra/code/dist/js/entity");
const RuntimeItemsMixin_1 = require("@mat3ra/code/dist/js/entity/mixins/RuntimeItemsMixin");
const flavorMixin_1 = require("./flavorMixin");
class Flavor extends entity_1.NamedDefaultableInMemoryEntity {
}
exports.default = Flavor;
// Apply mixins
(0, flavorMixin_1.flavorMixin)(Flavor.prototype);
(0, RuntimeItemsMixin_1.runtimeItemsMixin)(Flavor.prototype);
(0, flavorMixin_1.flavorStaticMixin)(Flavor);
