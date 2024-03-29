import { ContextProvider } from "@mat3ra/code/dist/js/context";

export class ExecutableContextProvider extends ContextProvider {
    constructor(config) {
        super({
            ...config,
            domain: "executable",
        });
    }
}
