"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.executableMixin = executableMixin;
function executableMixin(item) {
    // @ts-ignore
    const properties = {
        // get defaultFlavor() {
        //     return this.getFlavorByName();
        // },
        // get flavorsTree() {
        //     return this.prop("flavors", {}) as Record<string, any>;
        // },
        // get flavorsFromTree() {
        //     return Object.keys(this.flavorsTree).map((key) => {
        //         return (this.constructor as unknown as ExecutableStaticProperties).constructFlavor({
        //             ...this.flavorsTree[key],
        //             name: key,
        //         });
        //     });
        // },
        // getFlavorByName(name?: string) {
        //     return this.getEntityByName(this.flavors, "flavor", name || "") as BaseFlavor;
        // },
        // getFlavorByConfig(config?: { name: string }) {
        //     return config ? this.getFlavorByName(config.name) : this.defaultFlavor;
        // },
        // getFlavorsByApplicationVersion(version: string) {
        //     const filteredFlavors = this.flavors.filter((flavor) => {
        //         const supportedApplicationVersions = flavor.prop<string[]>(
        //             "supportedApplicationVersions",
        //         );
        //         return (
        //             !supportedApplicationVersions || supportedApplicationVersions.includes(version)
        //         );
        //     });
        //     return filteredFlavors;
        // },
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
