"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.materialsSetContextMixin = materialsSetContextMixin;
exports.MaterialsSetContextMixin = MaterialsSetContextMixin;
const utils_1 = require("@mat3ra/code/dist/js/entity/set/ordered/utils");
function materialsSetContextMixin(item) {
    const properties = {
        _materialsSet: undefined,
        get materialsSet() {
            return this._materialsSet;
        },
        initMaterialsSetContextMixin() {
            var _a;
            // @ts-ignore
            this._materialsSet =
                (_a = this.config.context) === null || _a === void 0 ? void 0 : _a.materialsSet;
        },
        sortMaterialsByIndexInSet(materials = []) {
            // DO NOT SORT IN PLACE AS IT CHANGES THE ORDER IN `this.materials` AND HAS SIDE EFFECTS (MaterialViewer).
            return materials.concat().sort((a, b) => {
                return (0, utils_1.compareEntitiesInOrderedSetForSorting)(
                    a,
                    b,
                    this.materialsSet._id,
                    false,
                );
            });
        },
    };
    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
}
function MaterialsSetContextMixin(superclass) {
    materialsSetContextMixin(superclass.prototype);
    return superclass;
}
