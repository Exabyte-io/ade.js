import { expect } from "chai";

import ContextProvider, { ContextProviderName } from "../../src/js/context/ContextProvider";

describe("ContextProvider", () => {
    const minimal = { name: ContextProviderName.KGridFormDataManager };
    const data = { a: "test" };

    it("can be created", () => {
        const provider = new ContextProvider(minimal);
        // eslint-disable-next-line no-unused-expressions
        expect(provider).to.exist;
    });

    it("sets and gets data", () => {
        const provider = new ContextProvider(minimal);
        provider.setData(data);
        expect(() => provider.getData()).to.throw("Not implemented.");
        provider.setIsEdited(true);
        expect(JSON.stringify(provider.getData())).to.equal(JSON.stringify(data));
        expect(() => provider.defaultData).to.throw("Not implemented.");
    });

    // transform, yieldData, yieldDataForRendering
});
