"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.executableMixin = executableMixin;
function executableMixin(item) {
    // @ts-ignore
    const properties = {
        toJSON(exclude = []) {
            const thisProto = Object.getPrototypeOf(this);
            const superProto = Object.getPrototypeOf(thisProto);
            const baseToJSON = superProto.toJSON;
            return baseToJSON.call(this, ["flavors"].concat(exclude));
        },
    };
    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
    return item;
}
