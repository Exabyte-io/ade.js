"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const entity_1 = require("@mat3ra/code/dist/js/entity");
const applicationMixin_1 = require("./applicationMixin");
const tree_1 = require("./tree");
class Application extends entity_1.NamedDefaultableInMemoryEntity {
    constructor(config) {
        const staticConfig = (0, tree_1.getApplicationConfig)(config);
        super({ ...staticConfig, ...config });
    }
}
exports.default = Application;
(0, applicationMixin_1.applicationMixin)(Application.prototype);
(0, applicationMixin_1.applicationStaticMixin)(Application);
