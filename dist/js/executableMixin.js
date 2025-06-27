"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.executableMixin = executableMixin;
function executableMixin(item) {
    // @ts-expect-error
    const properties = {
        get applicationId() {
            return this.prop("applicationId", "");
        },
        set applicationId(value) {
            this.setProp("applicationId", value);
        },
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
