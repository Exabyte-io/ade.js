"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.executableMixin = executableMixin;
function executableMixin(item) {
    // @ts-expect-error
    const properties = {
        get applicationId() {
            return this.prop("applicationId", []);
        },
        set applicationId(value) {
            this.setProp("applicationId", value);
        },
    };
    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
    return item;
}
