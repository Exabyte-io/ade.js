"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.executableMixin = executableMixin;
exports.executableStaticMixin = executableStaticMixin;
const flavor_1 = __importDefault(require("./flavor"));
function executableMixin(item) {
    // @ts-ignore
    const properties = {
        get defaultFlavor() {
            return this.getFlavorByName();
        },
        get flavorsTree() {
            return this.prop("flavors", {});
        },
        get flavors() {
            const tree = this.flavorsTree || {};
            return Object.keys(tree).map((key) => {
                return this.constructor.constructFlavor({
                    ...tree[key],
                    name: key,
                    executable: this,
                });
            });
        },
        get flavorsFromTree() {
            return Object.keys(this.flavorsTree).map((key) => {
                return this.constructor.constructFlavor({
                    ...this.flavorsTree[key],
                    name: key,
                });
            });
        },
        getFlavorByName(name) {
            return this.getEntityByName(this.flavors, "flavor", name || "");
        },
        getFlavorByConfig(config) {
            return config ? this.getFlavorByName(config.name) : this.defaultFlavor;
        },
        getFlavorsByApplicationVersion(version) {
            const filteredFlavors = this.flavors.filter((flavor) => {
                const supportedApplicationVersions = flavor.prop("supportedApplicationVersions");
                return (!supportedApplicationVersions || supportedApplicationVersions.includes(version));
            });
            return filteredFlavors;
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
function executableStaticMixin(Executable) {
    const properties = {
        constructFlavor(config) {
            if (this.constructCustomFlavor) {
                return this.constructCustomFlavor(config);
            }
            return new flavor_1.default(config);
        },
    };
    Object.defineProperties(Executable, Object.getOwnPropertyDescriptors(properties));
    return properties;
}
