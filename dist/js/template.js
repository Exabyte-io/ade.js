"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const entity_1 = require("@mat3ra/code/dist/js/entity");
const templateMixin_1 = require("./templateMixin");
class Template extends entity_1.NamedInMemoryEntity {
}
exports.default = Template;
// Apply mixins
(0, templateMixin_1.templateMixin)(Template.prototype);
(0, templateMixin_1.templateStaticMixin)(Template);
